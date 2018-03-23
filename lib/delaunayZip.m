function [F,V,C]=delaunayZip(F1,V1,F2,V2,inputStruct)

% function [Fn]=delaunayZip(F1,V1,F2,V2,inputStruct)
%-------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------

%%
[F,V,C]=joinElementSets({F1,F2},{V1,V2});
D=patchEdgeLengths(F,V);
[N,Vn,Nv]=patchNormal(F,V);

%%

defaultInputStruct.ind1=[];
defaultInputStruct.ind2=[];
defaultInputStruct.distLocal=2*max(D);
defaultInputStruct.startInd=[];
defaultInputStruct.plotOn=0;
[inputStruct]=structComplete(inputStruct,defaultInputStruct,0); %Complement provided with default if missing or empty

ind1=inputStruct.ind1;
ind2=inputStruct.ind2;
distLocal=inputStruct.distLocal;
startInd=inputStruct.startInd;
plotOn=inputStruct.plotOn;

%%

ind2=ind2+size(V1,1);
startInd(2)=startInd(2)+size(V1,1);

%%
%Create edges list
E=[ind1(1:end-1) ind1(2:end); ind2(1:end-1) ind2(2:end)];

%%

logicNotUsed=ismember((1:1:size(V,1))',E); %Logic to keep track of points that are used

%Initiate plot handles if plotting is on
if plotOn==1
    h1=[];
    h2=[];
    h3=[];    
end

%Initialize other parameters
Fn=[]; %Faces
Cn=[]; %"Colors"=step count for face group
c=1; %While loop counter variable

indGroup=startInd(:); %Initialize current group
numGroup=numel(indGroup); %Initialize current number of members of the current group

if plotOn==1
    cFigure; hold on;
    gpatch(F1,V1,'r','none',0.2);
    gpatch(F2,V2,'b','none',0.2);
    plotV(V(ind1,:),'r.-','LineWidth',1,'MarkerSize',10);
    plotV(V(ind2,:),'b.-','LineWidth',1,'MarkerSize',10);
    
    axisGeom;
    camlight headlight;
    colormap(gjet(250));
    colorbar;
    drawnow;
end

while 1
    
    if plotOn==1 %%Plot if plotting is on
        delete(h1);
        delete(h2);
        delete(h3);
    end
    
    if plotOn==1 %%Plot if plotting is on
        h3(2)=plotV(V(startInd,:),'y.','MarkerSize',60);
        if ~isempty(Fn)
            h3(3)=gpatch(Fn,V,'kw','k',0.5);
            h3(4)=patchNormPlot(Fn,V);
            h3(3)=plotV(V(logicNotUsed,:),'go','MarkerSize',20);
        end
        drawnow;
    end
    
    %Grow the current rebion
    while 1 %Loop to form local group (stuff attached to current group within a given distance)
        logicMember=any(ismember(E,indGroup),2); %Logic for all edges touching the current groupt
        E_sub=E(logicMember,:); %The subset of touching edges
        
        indGroup=unique([indGroup; E_sub(:)]); %Grow group with point indices in the edges that are touching
        [D]=minDist(V(indGroup,:),V(startInd,:));%Compute distance of group members to the start points
        logicKeep= (D<distLocal) & (logicNotUsed(indGroup));
        indGroup=indGroup(logicKeep); %Remove points that are too far
        if numGroup==numel(indGroup) %Compare current group size to previous step
            break %break while loop if the group is no longer growing
        end
        numGroup=numel(indGroup); %Get new current group size
    end
    
    if numel(indGroup)>2
        
        %Get curve start and end points
        logicMember=all(ismember(E,indGroup),2); %Logic for all edges touching the current groupt
        E_sub=E(logicMember,:); %The subset of touching edges
        
        [~,~,~,vCount]=cunique(E_sub); %Get vertex occurance counts
        indEndPoints=unique(E_sub(vCount==1));
        indEndPoints1=indEndPoints(ismember(indEndPoints,ind1));
        indEndPoints2=indEndPoints(ismember(indEndPoints,ind2));
        
        %Compose sub-curves
        logicSub1=all(ismember(E_sub,ind1),2);
        E_sub1=E_sub(logicSub1,:);
        [indListSub1]=edgeListToCurve(E_sub1);
        indListSub1=indListSub1(:);
        
        logicSub2=all(ismember(E_sub,ind2),2);
        E_sub2=E_sub(logicSub2,:);
        [indListSub2]=edgeListToCurve(E_sub2);
        indListSub2=indListSub2(:);
        
        %Create single closed curve
        if isempty(Fn) %First iteration
            [~,indClosest]=minDist(V(indListSub1(1),:),V([indListSub2(1) indListSub2(end)],:));
            if indClosest==2
                indListSub=[indListSub1;indListSub2;indListSub1(1)];
            else
                indListSub=[indListSub1;flipud(indListSub2);indListSub1(1)];
            end
        else
            indStart1=find(ismember(indListSub1,startInd));
            if indStart1~=1
                indListSub1=flipud(indListSub1);
            end
            indStart2=find(ismember(indListSub2,startInd));
            
            if indStart2~=1
                indListSub2=flipud(indListSub2);
            end
            
            indListSub=[indListSub1;flipud(indListSub2);indListSub1(1)];
            
        end
        indGroup=indListSub(1:end-1);
        
        if plotOn==1 %%Plot if plotting is on
            h1(1)=plotV(V(indEndPoints1,:),'r.','MarkerSize',50);
            h1(2)=plotV(V(indEndPoints2,:),'b.','MarkerSize',50);
            h1(3)=plotV(V(indListSub1,:),'r-','LineWidth',3);
            h1(4)=plotV(V(indListSub2,:),'b-','LineWidth',3);
            h1(5)=plotV(V(indListSub,:),'k-','LineWidth',2);
            drawnow;
        end
        
        %Rotate current point set
        V_now=V(indListSub(1:end-1),:); %Current closed curve coordinate set
        [R]=pointSetPrincipalDir(V_now); %Fit local coordinate system with 3rd direction pointing outward of local planar-ish region
        V_now_R=V_now*R; %Rotate coordinate set to prepare for 2D Delaunay based triangulation
        
        %Do 2D Delaunay triangulation
        DT_contraints=[(1:numel(indGroup))' ([2:numel(indGroup) 1])']; %Constraints are edges forming boundary
        DT = delaunayTriangulation(V_now_R(:,[1 2]),DT_contraints); % Initial triangulation
        f=DT.ConnectivityList; %Get faces set
        L = isInterior(DT); %Remove faces not inside region
        f=f(L,:); %Faces excluding external faces (only keep those inside constraint edges)
        f=indGroup(f); %Change indices to overall system
        
%         %Remove first and last triangles if there are more than 3 triangles, this
%         %will improve triangulation quality as these triangles were not fully
%         %embedded in the point set.
%         if size(f,1)>2
%             logicKeep=~(any(ismember(f,indEndPoints1(~ismember(indEndPoints1,Fn))),2) & any(ismember(f,indEndPoints2(~ismember(indEndPoints2,Fn))),2));
%             f=f(logicKeep,:);
%             
%             [IND_F,IND_V,IND_FF]=tesIND(f,V);
%             logicNeighbours=sum(IND_FF>1,2)>1;
%             logicNeighbours(any(ismember(f,startInd),2))=1;
%             f=f(logicNeighbours,:);            
%         end        
    else        
        f=unique([indGroup(:); startInd(:)])';       
    end

    if dot(mean(patchNormal(f,V),1),mean(Nv(indListSub(1:end-1),:),1))<0
        f=fliplr(f);        
    end

    if plotOn==1 %%Plot if plotting is on
        V_now_mean=mean(V_now); %Mean of coordinate set
        h2(1)=plotV(V_now_mean,'kx','MarkerSize',15);
        h2(2)=gpatch(f,V,'r','k',0.5);
        h2(3)=patchNormPlot(f,V);
        drawnow;
    end
    
    %Collect faces and color data
    Fn=[Fn;f];
    Cn=[Cn;c*ones(size(f,1),1)];
    
    %Get new start points
    indUsedNow=unique(Fn);

    %%
    E_Fn=sort(patchEdges(Fn,0),2);
    ind_E_Fn= reshape(sub2indn(size(V,1)*ones(1,2),E_Fn),size(Fn));
    
    Es=sort(E,2);
    ind_E= sub2indn(size(V,1)*ones(1,2),Es);
    
    logicEdgesUsed=ismember(ind_E,ind_E_Fn);
    
    E_sub=E(~logicEdgesUsed,:);
    indUnusedPoints=unique(E_sub);
    logicNotUsed=false(size(logicNotUsed));
    logicNotUsed(indUnusedPoints)=1;    
    
    if any(logicNotUsed)==0
        break
    end
    
    %Get curve start and end points
    [~,~,~,vCount]=cunique(E_sub); %Get vertex occurance counts
    indEndPoints=unique(E_sub(vCount==1));
         indEndPoints1=indEndPoints(ismember(indEndPoints,ind1));
        indEndPoints2=indEndPoints(ismember(indEndPoints,ind2));
        
    e=patchEdges(Fn,1);
    edgeEndPoints=e(any(ismember(e,indEndPoints1),2)&any(ismember(e,indEndPoints2),2),:);
    logicMemberOfLast=any(ismember(edgeEndPoints,f),2);
    edgeEndPoints=edgeEndPoints(logicMemberOfLast,:);
    if ~isempty(edgeEndPoints)
        startInd=edgeEndPoints(1,:);
    else
        startInd=[indListSub1(end) indListSub2(end)];
    end

        %%
    if plotOn==1 %%Plot if plotting is on
        h2(3)=plotV(V(logicNotUsed,:),'go','MarkerSize',20);
        h2(4)=plotV(V(startInd,:),'y.','MarkerSize',60);
    end
    
    c=c+1;

end

F=[F;Fn];
C=[C; max(C(:))+Cn];


