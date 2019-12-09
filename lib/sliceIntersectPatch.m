function [V_sliceContours]=sliceIntersectPatch(F,V,N_slices,V0_slices,plotOptions)

markerSize=25; 

maxEpsV=max(eps(V(:)));

if isempty(plotOptions)
    plotOn=0;
else
    plotOn=1;
    fontSize=20;
    
    hf1=cFigure;
    title('slicePatch','FontSize',fontSize);    
    hold on;
    
    hp= patch('Faces',F,'Vertices',V,'FaceColor',plotOptions.faceColorPatch,'EdgeColor','k','FaceAlpha',plotOptions.faceAlphaPatch);
    
    camlight('headlight'); lighting flat;
    axisGeom;
    drawnow;
    
    if plotOptions.plotPlanes
        axisLim=axis;
        %Creating slice plane for plotting
        F_slicePlane=[1 2 3 4];
    end
%     if isfield(plotOptions,'saveFrames')
%         saveOn=1;
%     else
%         saveOn=0;
%     end
    
end

numSlices=size(V0_slices,1);

if plotOn==1;
    h1=[]; h2=[]; h3=[]; h4=[];
end
V_sliceContours=cell(numSlices,1);

splitGroupDetect=0;
for q=1:1:numSlices;
    
    %% Create rotation matrix
    sliceNormal=N_slices(q,:); %Current slice normal
    e3=sliceNormal; e3=vecnormalize(e3);
    dcmDone=0; runInc=1;
    while dcmDone==0 %While loop to create valid orthonormal system based on slice normal
        switch runInc;
            %Try coordinate system base vectors first
            case 1
                e1=[1 0 0];
            case 2
                e1=[0 1 0];
            case 3
                e1=[0 0 1];
            otherwise %a random vector
                e1=randn(1,3); e1=vecnormalize(e1);
        end
        
        e2=cross(e3,e1);
        if ~all(e2<maxEpsV) && ~all(isnan(e2))%If the random vector is not colinear with the normal
            dcmDone=1;
        end
        e2=vecnormalize(e2);
        
        if runInc==10
            fdsaf
        end
        runInc=runInc+1;
    end
    e1=cross(e2,e3); e1=vecnormalize(e1);
    DCM=[e1; e2; e3]; %The rotation or direction cosine matrix
    
    %% Calculate intersections
    
    %Rotating to Z-slice problem
    V0_slices_rot=V0_slices/DCM; %rotated cut points
    V_rot=V/DCM; %rotated patch coordinates
    
    Z_cut=V0_slices_rot(q,3); %Slice level
    
    % Plotting
    if plotOn==1 && plotOptions.plotPlanes;
        if plotOptions.cleanUp==1
            delete(h1);
        end
        V_slicePlane=[axisLim(1) axisLim(3) Z_cut;...
            axisLim(2) axisLim(3) Z_cut;...
            axisLim(2) axisLim(4) Z_cut;...
            axisLim(1) axisLim(4) Z_cut];
        h1= patch('Faces',F_slicePlane,'Vertices',V_slicePlane*DCM,'FaceColor','b','EdgeColor','w','FaceAlpha',0.25);
    end
        
    Z=V_rot(:,3); %Z coordinates of all rotated vertices
    
    %Slightly offset coordinates that coincide with cut level
    %TO DO: Fix this in the future 
    Z_offset=eps(Z_cut);
    V_rot(Z==Z_cut,3)=V_rot(Z==Z_cut,3)+Z_offset;    
    Z=V_rot(:,3); %Overwrite Z coordinates of all rotated vertices        
            
    ZF=Z(F); %Z coordinates of face vertices
   
    logicBelow=(ZF<=Z_cut);
    logicBelowSumTwo=sum(logicBelow,2)==2;
    
    logicAbove=(ZF>=Z_cut);
    logicAboveSumTwo=sum(logicAbove,2)==2;
    
    logicTriValid=(logicBelowSumTwo | logicAboveSumTwo) & (any(logicBelow,2) & any(logicAbove,2)); %Logic for faces have some vertices on or below the cutting plane and above
    
    if any(logicTriValid) %if none are cut we skip
        
        %Selecting triangles that are cut
        Fcut=F(logicTriValid,:);
        logicBelowCut=logicBelow(logicTriValid,:);
        logicAboveCut=logicAbove(logicTriValid,:);
        
        triGroups=tesgroup(Fcut); %Grouping to separate triangle sets
        numTriGroups=size(triGroups,2);
        
        if plotOn==1;
            if isempty(plotOptions.faceColorSubPatch)
                patchColors=jet(numTriGroups+1);
            else
                patchColors=repmat(plotOptions.faceColorSubPatch,numTriGroups,1);
            end
        end
        
        curveGroupCount=1; %N.B. this may differ from size(triGroups,2)
        for groupId=1:1:size(triGroups,2);
            
            %Selecting current group
            logicGroup=triGroups(:,groupId);
            Fsub=Fcut(logicGroup,:);
            
            logicAboveCutSub=logicAboveCut(logicGroup);
            logicBelowCutSub=logicBelowCut(logicGroup);
            
            %Plotting
            if plotOn==1;
                if plotOptions.cleanUp==1
                    try delete(h2); delete(h3); end
                end
                h2= patch('Faces',Fsub,'Vertices',V,'FaceColor',patchColors(groupId,:),'EdgeColor','w','FaceAlpha',1);
            end
            
            %Get edges
            E_1=[Fsub(:,1) Fsub(:,2)];
            E_2=[Fsub(:,2) Fsub(:,3)];
            E_3=[Fsub(:,3) Fsub(:,1)];
            E=sort([E_1; E_2; E_3],2);            
           
            indEdges=sub2ind(max(F(:))*ones(1,2),E(:,1),E(:,2)); %Virtual edge indices
            indEdgesPerFace=reshape(indEdges,size(Fsub,1),3);
            
            [~,indUni,~]=unique(indEdges);
            E=E(indUni,:); %The unique edge set
            indEdges=indEdges(indUni); %The unique virtual index set
                        
            %Remove invalid edges which produce no intersection
            logicEdgesInvalid=all(Z(E)<Z_cut,2) | all(Z(E)>Z_cut,2); %All above or below
            E=E(~logicEdgesInvalid,:); %The valid edge set
            indEdges=indEdges(~logicEdgesInvalid); %The valid virtual index set
            
            %Calculate intersection points
            V0=(ones(size(E,1),1)*V0_slices_rot(q,:)); %Arbitrary point on slice plane
            P1=V_rot(E(:,1),:); %first point
            P2=V_rot(E(:,2),:); %second point
            D1=V0-P1; D2=P2-P1; %Difference measures
            dotND1=D1(:,3); %Simplified dot product due to [0 0 Zc]
            dotND2=D2(:,3);
            R1=(dotND1./dotND2);
            Vi=P1+(R1*ones(1,size(Fsub,2))).*D2;
            
%             %Fix intersection for colinear edges
%             L_colinear=(dotND2==0);            
%             Xcor=mean([P1(L_colinear,1) P2(L_colinear,1)],2);
%             Ycor=mean([P1(L_colinear,2) P2(L_colinear,2)],2);
%             Zcor=mean([P1(L_colinear,3) P2(L_colinear,3)],2);
%             Vi(L_colinear,:)=[Xcor Ycor Zcor];           
            Vi(:,3)=Z_cut; 
            
            if plotOn==1
                if plotOptions.cleanUp==1
                    try delete(h1); end
                    try delete(h2); end
                    try delete(h3); end
                    try delete(h4); end
                end
                Vp=Vi*DCM;
                h1=plot3(Vp(:,1),Vp(:,2),Vp(:,3),'b.','MarkerSize',markerSize);
            end
            
            %% START WHILE LOOP
            
            done=0;
            qw=1;
            indFace_now=1;
            logicFacesDone=false(size(Fsub,1),1);
            logicEdgesDone=false(size(E,1),1);
            V_path=nan(size(Fsub,1),3);
            
            trySplit=0;
            clc;
            while done==0
                
                %Setting current parameters
%                 F_now=Fsub(indFace_now,:); %Current face

                indEdgesInd_now=indEdgesPerFace(indFace_now,:);
                
                %                 if ~isempty(indEdgesInd_now);
                
                L_toDo=ismember(indEdgesInd_now,indEdges(~logicEdgesDone));               
                        
                indEdgeInd_now=indEdgesInd_now(L_toDo); %current edge index
                if ~isempty(indEdgeInd_now)
                    if qw==1
                        indEdgeInd_now=indEdgeInd_now(1); %first time the first one is picked
                    end
                    if numel(indEdgeInd_now)>1 %Cross roads?
                        if qw==2
                            indEdgeInd_now=indEdgeInd_now(1); %Pick first again
                            indEdge_now=find(indEdges==indEdgeInd_now);
                        else %Cross roads after some patch elements have been established now you could choose based on path history
                            indEdgeInd_now=indEdgeInd_now(1); %Pick first again for now
                            indEdge_now=find(indEdges==indEdgeInd_now);
                            disp('Check for crossroads');
                        end
                    else %Normal case, just 1 option
                        indEdge_now=find(indEdges==indEdgeInd_now);
                    end
                    
                    distPoints=abs(V_path-(ones(size(V_path,1),1)*Vi(indEdge_now,:)));
                    logicDistPoints=all(distPoints<maxEpsV,2);
                    
                    if all(logicDistPoints==0) %If point is truely new
                        V_path(qw,:)=Vi(indEdge_now,:); %add current point to patch
                        V_sliceContours{q}{curveGroupCount}(qw,:)=Vi(indEdge_now,:)*DCM; %add current point to slice contours
                        qw=qw+1; %Increment counter
                    end
                    
                    %Plotting
                    if plotOn==1 && plotOptions.animate==1;
                        if plotOptions.cleanUp==1
                            try delete(h4); end
                        end
                        Vp=V_path*DCM;
                        h4=plot3(Vp(:,1),Vp(:,2),Vp(:,3),'g.-','LineWidth',10,'MarkerSize',markerSize);
                        drawnow;
                    end
                    
                    %Remove treated values
                    logicEdgesDone(indEdge_now)=1;
                    logicFacesDone(indFace_now)=1;
                    
                    %Find other (not current) face containing this edge
                    logicFaceNext=any(ismember(indEdgesPerFace,indEdges(indEdge_now)),2)&~logicFacesDone;                    
                                    
                    if any(logicFaceNext)
                        indFace_now=find(logicFaceNext);
                    elseif ~all(logicEdgesDone)
                        indFace_now=find(~logicFacesDone,1);
                        splitGroupDetect=1;
                        disp('Curve path splits face group. Creating additional curve set (1)');                     
                    end
                                        
                elseif ~all(logicEdgesDone)
                    indFace_now=find(~logicFacesDone,1);
                    splitGroupDetect=1;
                    disp('Curve path splits face group. Creating additional curve set (2)');
                    
                    if trySplit==1; %Already tried that, move on
                        done=1;
                        disp('STOPPED WHILE LOOP, SPLITTING UNSUCCESSFUL');
                    end
                    trySplit=1;
                end
                
                
                if splitGroupDetect==1;
                    qw=1; %Reset increment to 1
                    V_path=nan(nnz(~logicEdgesDone),3); %Reset path
                    splitGroupDetect=0;
                    curveGroupCount=curveGroupCount+1;
                end
                
                if all(logicEdgesDone)
                    done=1;
                end
                
                %                 else
                %                     done=1;
                %                 end
            end
            %%
            
            if plotOn==1;
                if plotOptions.cleanUp==1
                    try delete(h3); end
                    try delete(h4); end
                end
                
                if curveGroupCount~=groupId
                    Vp=[];
                    for qp=1:1:numel(V_sliceContours{q})
                        Vp=[Vp; nan(1,3); V_sliceContours{q}{qp}];
                    end
                else
                    Vp=V_path*DCM;
                end
                
                h3=plot3(Vp(:,1),Vp(:,2),Vp(:,3),'g.-','LineWidth',5,'MarkerSize',markerSize); drawnow;
            end
            if plotOn==1
                pause(plotOptions.pauseTime)
            end
            curveGroupCount=curveGroupCount+1;
        end
    end
%     if plotOn==1 && saveOn==1;
%         saveName=[plotOptions.saveFrames(1:end-4),'_',num2str(q),plotOptions.saveFrames(end-3:end)];
%         axis off;
%         export_fig(saveName,'-r100');
%         axis on;
%     end
end

if plotOn==1;
    if plotOptions.cleanUp==1
        try delete(h1); end
        try delete(h2); end
        try delete(h3); end
        try delete(h4); end
        close(hf1);
    end
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
