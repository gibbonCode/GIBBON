function [Vcs]=uiContourSegment(varargin)

% function [Vcs]=uiContourSegment(M,cPar,saveName)
% ------------------------------------------------------------------------
% Uppon running this function a basic user interface segmentation program
% is started. Users can step through the slices of the image volume M and
% define contours for each slice segmenting features of interest.
%
% By pressing H during operation the following help information appears on
% screen: 
%
%                     'Up              -> Increase contour level',...
%                     'Down            -> Decrease contour level',...
%                     'Left click      -> Select closest contour to keep',...
%                     'Right click     -> Delete closest contour',...
%                     'C               -> Define crop window for contours',...
%                     'D               -> Manually draw points as a contour, left click to draw, right to finish',...
%                     'V               -> Add selected contours from the previous slice for use in current slice',...
%                     'O               -> Reorder contour points. Warning this is a slow process work in ordered e.g. clockwise fashion to avoid ordering',...
%                     'R               -> Remove all selected/drawn contours',...
%                     'S               -> Create split contour, i.e. keep current contour segment but add aditional for this slice',...
%                     'Space bar       -> Done with current slice, proceed to next',...
%                     'H               -> Display help',...                    
%
%
% Added basic documentation text: 23/10/2013
% TO DO: Add more documentation. Turn into an actual GUI
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/10/23
% 2016/02/25 Fixed bug in relation to padding when combined with a user
% specified slice range. 
% ------------------------------------------------------------------------

%%

switch nargin
    case 3
        M=varargin{1};
        cPar=varargin{2};
        saveName=varargin{3};        
        guideContourCell={};
    case 4        
        M=varargin{1};
        cPar=varargin{2};
        saveName=varargin{3};
        guideContourCell=varargin{4};
    otherwise
        error('Wrong number of input arguments');
end

%%
siz=size(M);

%Cope with 2D image sets
if ismatrix(M)
    siz(3)=1;
end

%% PLOT SETTINGS

figStruct.Name='GIBBON'; %Figure name
figStruct.Color='k'; %Figure background color
figStruct.ColorDef='black'; %Setting colordefinitions to black
figStruct.ScreenOffset=0; %Setting spacing of figure with respect to screen edges
fontSize=10;
mark_siz1=2;
mark_siz2=50;
lineWidth1=1;
lineWidth2=2;
F_alpha1=1;
cMap=gray(250);

%% GET CONTROL PARAMETERS

%Only contours with more points are considered
if isfield(cPar,'minContourSize')
    minContourSize=cPar.minContourSize;
else
    minContourSize=150; %DEFAULT
end

%Degree of smoothening of csaps function
if isfield(cPar,'smoothFactor')
    smoothFactor=cPar.smoothFactor;
else
    smoothFactor=0.25; %DEFAULT
end

%Reduction factor for contour smoothening
if isfield(cPar,'pointReductionFactor')
    pointReductionFactor=cPar.pointReductionFactor;
else
    pointReductionFactor=5; %DEFAULT
end

%Ones describe image data regions of interest
if isfield(cPar,'logicBackGround')
    logicBackGround=cPar.logicBackGround;
    if isempty(logicBackGround)
        logicBackGround=true(siz); %DEFAULT
    end
else
    logicBackGround=true(siz); %DEFAULT
end

%Check padOn option
if isfield(cPar,'padOn')
    padOn=cPar.padOn;
else
    padOn=1; %DEFAULT
end

%Check recoverOn options
if isfield(cPar,'recoverOn')
    recoverOn=cPar.recoverOn;
else
    recoverOn=0; %DEFAULT
end

%Get voxel size
if isfield(cPar,'v')
    v=cPar.v;
else
    v=[1 1 1]; %DEFAULT
end

%Get slice range
if isfield(cPar,'sliceRange')
    sliceRangeUser=cPar.sliceRange;
else
    ind=find(logicBackGround);
    [~,~,k]=ind2sub(size(logicBackGround),ind);
    sliceRangeUser=sort(unique(k(:))); %DEFAULT
end

%Distance thresholds for curve spacing, higher spacings are considered gaps for instance
distThresh=2*max(v(1:2));

%% Pad arrays to cope with open ended contours
if padOn==1
    
    M_p=zeros(siz(1)+2,siz(2)+2,siz(3)); %Padded in row and column direction
    M_p(2:end-1,2:end-1,:)=M; %Assign image data in the middle
    M=M_p; %Override M
    siz=size(M); %Overide size
    
    %Fix background logic
    logicBackGround_p=false(size(M));
    logicBackGround_p(2:end-1,2:end-1,:)=logicBackGround;
    logicBackGround=logicBackGround_p; %Overide background logic    
    
end

%% SET-UP IMAGE AND IMAGE COORDINATES

%Normalizing the image data
M=dataNorm(M,1);

%Setting background to zero
M(~logicBackGround)=0;

%Get image coordinates
[J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));

%Convert to cartesian coordinates using voxels size if provided
[X,Y,Z]=im2cart(I,J,K,v);


%% Creating slice range and logic for slice plotting

%Based on background logic if provided
IND_BG = find(logicBackGround);
[~,~,Kbg] = ind2sub(siz,IND_BG);
sliceRange=unique(Kbg(:))'; %Slice range based on background logic
logicSliceRange=ismember(sliceRange,sliceRangeUser);
sliceRange=sliceRange(logicSliceRange); %Slices to analyse based on input and logic
sliceMinInd=min(sliceRange(:)); %sliceMaxInd=max(sliceRange(:));
L_slice=any(logicBackGround,3); %logic for slice plotting

%%

sliceMinInd

%Plotting first slice
L=false(siz);
L(:,:,sliceMinInd)=L_slice;
IND_slice=find(L);
[Fs,Vs,C_data]=ind2patch(IND_slice,M,'sk');
[Vs(:,1),Vs(:,2),Vs(:,3)]=im2cart(Vs(:,2),Vs(:,1),Vs(:,3),v);
Vs(:,3)=0; %The z-level for this function will be zero for all plotting

%Plotting image
hf1=cFigure(figStruct);
figNameString=' --- SEGMENTING CONTOURS --- ';
titleStringAdd=' ';
set(hf1,'Name',figNameString);

xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hp1=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none','CData',C_data,'FaceColor','flat','FaceAlpha',F_alpha1);
axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize);
colormap(cMap); colorbar; caxis([0 1]);
drawnow;

fdsadfas

setDefaultPointer; %Set default pointer
    
minC=0;
maxC=1;

Vcs=cell(1,siz(3)); 
Vcs{sliceMinInd}{1}=[]; %Result cell array
hc=[]; hvcs=[]; hvcs_guide=[]; %Initialise graphic handles
stepSize_tc=0.025; %Contour level step increment size
n_tc=round(1/stepSize_tc);
Tc=linspace(0,1,n_tc); %The contour level range
ic=round(numel(Tc)/6); %The initial contour index
z_offset=0.01; %Arbitrary offset ensuring curves are above slice during plotting

V=[];
% ic_old=ic; 

for qSlice=sliceRange
    qSlice
    
    %Check for recovery file
    recoveryFileName=['temp_recover_uiContourSegment_',num2str(qSlice),'.mat'];
    if recoverOn==1 %Initialize with recovered metrics
        if exist(recoveryFileName,'file')==2 %Initialize based on recovery file
            load(recoveryFileName);
            Vcs{qSlice}=recoveryStruct.Vcs;
            splitInd=recoveryStruct.splitInd;
            V=recoveryStruct.V;
            C=recoveryStruct.C;
            ic=recoveryStruct.ic;
            ic_old=recoveryStruct.ic_old;
            caxis([recoveryStruct.minC recoveryStruct.maxC]); %Set axis limits
        else
            Vcs{qSlice}{1}=[]; %Initialise first contour as empty for each slice
            splitInd=1; %Initialise to 1 for each slice
            ic_old=NaN;
        end
    else %Initialize normally
        Vcs{qSlice}{1}=[]; %Initialise first contour as empty for each slice
        splitInd=1; %Initialise to 1 for each slice        
        ic_old=NaN;
    end    
    
    %Plot new slices    
    %Redefine CData
    L=false(siz); 
    L(:,:,qSlice)=L_slice;
    [~,~,C_data]=ind2patch(L,M,'sk');        
    set(hp1,'CData',C_data); %Set FaceVertexCData
    zs=(qSlice-0.5).*v(3); %z-level for contourslice
    
    %Plot guide contours
    if ~isempty(guideContourCell)
        %Remove guide contour plot and handle
        if ~isempty(hvcs_guide)
            delete(hvcs_guide); hvcs_guide=[];
        end
        
        %Plot guide contours
        for q_guide=1:1:numel(guideContourCell)
            Vcs_guide=guideContourCell{q_guide};
            if ~isempty(Vcs_guide{qSlice})
                for qc=1:numel(Vcs_guide{qSlice})
                    if ~isempty(Vcs_guide{qSlice}{qc})
                        Vp=Vcs_guide{qSlice}{qc};
                        hvcsq_guide=plot3(Vp(:,1),Vp(:,2),z_offset.*ones(size(Vp(:,1))),'g-','lineWidth',lineWidth2); hold on;
                        hvcs_guide=[hvcs_guide hvcsq_guide];
                    end
                end                
            end
        end
        drawnow;                
    end
    
    done=0;
    while done==0; %While loop for current slice analysis
        
        if ic~=ic_old %Update contour only if level is changed
            
            [C,C_siz]=contour_group(X,Y,Z,M,[],[],zs,[Tc(ic) Tc(ic)]);
            
            C=C(C_siz>minContourSize); %Discarding contours that are too smalls
  
            %Smooth and resample contours
            if pointReductionFactor>1 || smoothFactor<1
                for q_smooth=1:numel(C)
                    
                    V=C{q_smooth};
                    nV=size(V,1);
                    %n=size(V,1);
                    
                    %Resampling
                    if pointReductionFactor>1
                        D=pathLength(V); %Compute distance metric used for parametric representation
                        n=round((max(D(:))./min(v))/pointReductionFactor);
                        %Limiting n to minContourSize
                        if n<minContourSize;
                            n=minContourSize; %use at least minContourSize points
                        end
                        V=V(1:end-1,:); %Trim off last point
                        [V]=evenlySampleCurve(V,n,'pchip',1);
                    else
                        V=V(1:end-1,:); %Trim off last point
                    end
           
                    %Smoothening
                    if smoothFactor<1
                        [V]=evenlySampleCurve(V,nV,smoothFactor,1);
                    end
                    C{q_smooth}=V;
                end
            end
            ic_old=ic;
        end
                
        %Plotting results for this slice        
        if ~isempty(hvcs)
            delete(hvcs); hvcs=[];
        end

        if ~isempty(Vcs{qSlice})
            hvcs=[];
            for qc=1:numel(Vcs{qSlice})
                if ~isempty(Vcs{qSlice}{qc})
                    Vp=Vcs{qSlice}{qc};
                    hvcsq=plot3(Vp(:,1),Vp(:,2),z_offset.*ones(size(Vp(:,1))),'b.-','MarkerSize',mark_siz1,'lineWidth',lineWidth2); hold on; 
                    hvcs=[hvcs hvcsq];
                end
            end
            axis equal; view(2); axis tight;  zlim([-1 1]); set(gca,'FontSize',fontSize);            
        end
        drawnow;
        
        %Plotting contours
        
        if ~isempty(hc)
            for ih=1:1:numel(hc)
                delete(hc(ih));
            end
            hc=[];
        end
        contourColorMap=autumn(numel(C));
        for ig=1:1:numel(C)
            V=C{ig};
            h=plot3(V(:,1),V(:,2),z_offset.*ones(size(V(:,3))),'k-','lineWidth',lineWidth1);
            set(h,'Color',contourColorMap(ig,:));
            hc=[hc; h];
        end
        axis equal; view(2); axis tight;  zlim([-1 1]); set(gca,'FontSize',fontSize);
        drawnow;
        title(['Contour slice:',num2str(qSlice),' Tc=',num2str(Tc(ic)),titleStringAdd]);
        
        %GET GINPUT
        [xc,yc,b]=qginput(1);

        switch b %Switch structure for various inputs
            case 106 %J key 
                %Change contour increment size
                prompt = {'Change contour increment size e.g. 0.05 (0-1):'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'0.05'};
                IC_answer = inputdlg(prompt,dlg_title,num_lines,def);
                if ~isempty(IC_answer)
                    stepSize_tc=str2double(IC_answer{1});                   
                end
                tcc=Tc(ic); %The current level                
                n_tc=round(1/stepSize_tc); %Number of steps
                Tc=linspace(0,1,n_tc); %The contour level range
                [~,ic]=min(abs(Tc-tcc)); %Adjusted index ic to keep current level 
            case 109 %M key
                prompt = {'Enter caxis min level:','Enter caxis max level:'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'0','1'};
                CAXIS_answer = inputdlg(prompt,dlg_title,num_lines,def);
                if ~isempty(CAXIS_answer)
                    minC=str2double(CAXIS_answer{1});
                    maxC=str2double(CAXIS_answer{2});
                    caxis([minC maxC]);
                end
            case 108 %L key
                prompt = {'Jump to contour level (0-1):'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'0.5'};
                IC_answer = inputdlg(prompt,dlg_title,num_lines,def);
                if ~isempty(IC_answer)
                    tcc=str2double(IC_answer{1});
                    [~,ic]=min(abs(Tc-tcc));
                end
            case 30 % Up key
                titleStringAdd=' Increased contour level';
                ic=((ic+1).*(ic+1<numel(Tc)))+numel(Tc).*(ic+1>=numel(Tc));
            case 31 %Down key
                titleStringAdd=' Decreased contour level';
                ic=((ic-1).*(ic-1>0))+(ic-1<1);
            case 99 % The C key -> Crop contours
                title('Cropping: Press and hold to drag and define a cropping window');
                
                %Draw cropping window
                
                %OLD METHOD
                %                 [xcc,ycc]=qginput(1,specialPointerShape('ulc'));
                %                 hcc1=plot3(xcc,ycc,z_offset,'r+','MarkerSize',mark_siz2);
                %                 [ic1,jc1,kc1]=cart2im(xcc,ycc,0,v);
                %                 [xcc,ycc]=qginput(1,specialPointerShape('lrc'));
                %                 hcc2=plot3(xcc,ycc,z_offset,'r+','MarkerSize',mark_siz2);
                %                 [ic2,jc2,kc2]=cart2im(xcc,ycc,0,v);
                
                %Using rbbox to draw a cropping window
                [mousePointerType]=specialPointerShape('cut');
                set(gcf,'Pointer','custom','PointerShapeCData',mousePointerType.PointerShapeCData,'PointerShapeHotSpot',mousePointerType.PointerShapeHotSpot); %User specified mousePointerType
                
                waitforbuttonpress;
                point1 = get(gca,'CurrentPoint');    % button down detected                
                rbbox;                 % return figure units
                drawnow; %pause(1e-6); % Required to avoid point1=point2 errors          
                point2 = get(gca,'CurrentPoint');    % button up detected                
                point1 = point1(1,1:2);            % extract x and y
                point2 = point2(1,1:2);
                p1 = [min(point1(1,1),point2(1,1)) max(point1(1,2),point2(1,2))];
                p2 = [max(point1(1,1),point2(1,1)) min(point1(1,2),point2(1,2))];
                [ic1,jc1,~]=cart2im(p1(1),p1(2),0,v);
                [ic2,jc2,~]=cart2im(p2(1),p2(2),0,v);
                hcc1=plot3(p1(1),p1(2),z_offset,'r+','MarkerSize',mark_siz2);
                hcc2=plot3(p2(1),p2(2),z_offset,'r+','MarkerSize',mark_siz2);
                drawnow;
                
                mc=false(siz(1:2));
                mc(round(ic2):round(ic1),round(jc1):round(jc2))=1;
                cropInd=find(mc);
                
                %Remove contour elements within crop window and regroup
                num_C=numel(C);
                L_delete=false(1,num_C);
                for ig=1:1:num_C
                    V=C{ig};
                    [icg,jcg,~]=cart2im(V(:,1),V(:,2),zeros(size(V(:,1))),v);
                    icg=round(icg);
                    jcg=round(jcg);
                    icg(icg<1)=1;
                    jcg(jcg<1)=1;
                    icg(icg>siz(1))=siz(1);
                    jcg(jcg>siz(2))=siz(2);
                    
                    vertexInd= sub2ind(siz(1:2),icg,jcg);
                    L_crop=~ismember(vertexInd,cropInd);
                    
                    if any(~L_crop(:))
                        if all(L_crop(:))
                            L_delete(ig)=1; %Whole contour will be removed
                        else %Contour is cropped
                            %Split curve in groups
                            groupIndices=cumsum(diff([0; L_crop]).*L_crop).*L_crop;
                            groupIndices(groupIndices==0)=NaN;
                            
                            D_startEnd=sqrt(sum((V(1,:)-V(end,:)).^2));
                            if all(~ismember([1 size(V,1)],cropInd)) && D_startEnd<distThresh% If the start and end are not cropped reunite start/end groups
                                groupIndices(groupIndices==1)=max(groupIndices(:));
                            end
                            
                            groupFirstIter=1;
                            for q_group=nanmin(groupIndices(:)):1:nanmax(groupIndices(:))
                                
                                Vg=V(groupIndices==q_group,:); %Current curve
                                
                                %Reorder so within curve gaps are removed
                                Dx=diff(Vg(:,1)); Dy=diff(Vg(:,2));
                                D=sqrt(Dx.^2+Dy.^2);
                                LD=D>distThresh;
                                if any(LD); %reorder
                                    indGap=find(LD,1)+1;
                                    Vg=[Vg(indGap:end,:); Vg(1:indGap-1,:)];
                                end
                                
                                if groupFirstIter==1;
                                    C{ig}=Vg;
                                    groupFirstIter=0;
                                else
                                    C{end+1}=Vg;
                                end
                            end
                        end
                    end
                end
                
                indDelete=find(L_delete);
                L_delete=false(1,numel(C));
                L_delete(indDelete)=1;
                
                C=C(~L_delete);
                
                %Get new group sizes
                group_sizes = cell2mat(cellfun(@(x) numel(x), C,'UniformOutput',0)');
                %Sort C array according to group size
                [C_siz,ind_sort]=sort(group_sizes);
                C=C(ind_sort);
                %Remove zero-size groups
                C=C(C_siz>0); C_siz=C_siz(C_siz>0);
                delete(hcc1); delete(hcc2);
                setDefaultPointer; %Set default pointer
                titleStringAdd=' Cropped contours';
            case 1 %Right click key -> Select group
                titleStringAdd=' Added selected contour to current slice';
                
                D=NaN(1,numel(C));
                for ig=1:1:numel(C)
                    V=C{ig};
                    D(ig)=min(hypot((V(:,1)-xc),(V(:,2)-yc)));
                end
                [~,indMin] = min(D);
                
                %Get the closest contour to the click
                V_add=C{indMin}; 
                
                %Remove contour from list
                Lg=true(1,numel(C)); Lg(indMin)=0;
                C=C(Lg);
                
                if ~isempty(Vcs{qSlice}{splitInd})
                    %Reorder selection before adding to maintain overall
                    %contour point order.
                    
                    %See if start or end is closest to current contour end
                    V_now=Vcs{qSlice}{splitInd}; %The contour so far
                    
                    V_now_start_end=V_now([1 size(V_now,1)],:);
                    V_add_start_end=V_add([1 size(V_add,1)],:);
                    try
                        D=dist(V_now_start_end,V_add_start_end');
                    catch
                        D=distND(V_now_start_end,V_add_start_end);
                    end
                    [~,indMin_lin]=min(D(:));
                    
                    switch indMin_lin
                        case 1 %1 1 -> starts are close, add in front and flip
                            V_add=flipud(V_add);
                            Vcs{qSlice}{splitInd}=[V_add; V_now];
                        case 2 %2 1 -> start is close to end, just add behind
                            Vcs{qSlice}{splitInd}=[V_now; V_add];
                        case 3 %1 2 -> end is close to start, just add in fron
                            Vcs{qSlice}{splitInd}=[V_add; V_now];
                        case 4 %2 2 -> ends are clsoe, add behind and flip
                            V_add=flipud(V_add);
                            Vcs{qSlice}{splitInd}=[V_now; V_add];
                    end
                else
                    Vcs{qSlice}{splitInd}=V_add; %Add selected group to Vcs{is}{splitInd}
                end
            case 111 %The O key -> Fix contour point order SLOW PROCESS
                title('Reordering points, may be slow, please wait...');
                [Vcs{qSlice}{splitInd},~]=curvePathOrderFix(Vcs{qSlice}{splitInd});
                titleStringAdd=' Reordered points';
            case 118 %The V key -> Add contour from previous slice
                titleStringAdd=' Added contour from previous slice';
                if qSlice>1
                    if ~isempty(Vcs{qSlice-1})                                                
                        for qc=1:numel(Vcs{qSlice-1})
                            if ~isempty(Vcs{qSlice-1}{qc})
                                V_add=Vcs{qSlice-1}{qc};
                                V_add(:,3)=zs; %Fix z coordinate for current slice
                                C{end+1}=V_add;
                            end
                        end
                    end
                end                                
            case 3 %Right click -> Delete nearest contours
                titleStringAdd=' Deleted nearest contour';
                D=NaN(1,numel(C));
                for ig=1:1:numel(C)
                    V=C{ig};
                    D(ig)=min(hypot((V(:,1)-xc),(V(:,2)-yc)));
                end
                [~,indMin] = min(D);
                Lg=true(1,numel(C)); Lg(indMin)=0;
                C=C(Lg);
            case 114 %The R key -> Reset/clear (remove all) selected contours
                titleStringAdd=' Resetted contours';
                Vcs{qSlice}{splitInd}=[];
            case 115 %The S key -> Split contour into cells
                titleStringAdd=' Splitting contour adding additional contour entry for this slice';
                splitInd=splitInd+1;
                Vcs{qSlice}{splitInd}=[];
            case 100 %The D key -> Draw points as seperate contour
                title('Draw manual points, left click to draw, right to finish');
                [mousePointerType]=specialPointerShape('pen');
%                 set(gcf,'Pointer','custom','PointerShapeCData',mousePointerType.PointerShapeCData,'PointerShapeHotSpot',mousePointerType.PointerShapeHotSpot); %User specified mousePointerType
                
                hpd=[];
                drawDone=0;
                Vd=[];
                while drawDone==0
                    [Xd,Yd,bd]=qginput(1,mousePointerType);
                    switch bd
                        case 1
                            Vd=[Vd; Xd Yd zs];
                            hpdi=plot3(Vd(:,1),Vd(:,2),z_offset*ones(size(Vd,1),1),'g.-','MarkerSize',mark_siz1,'lineWidth',lineWidth1);
                            hpd=[hpd hpdi]; %Collect graphic handle
                        case 3
                            drawDone=1;
                    end
                end
                delete(hpd);
                C{end+1}=Vd;
                setDefaultPointer;
                titleStringAdd=' Manual contour created';
            case 32 % Space bar to go to next slice
                done=1;                
            case 104 %The H key -> Display help
                msgText={'INPUT OPTIONS:',...
                    '------------------------------------------------------------------------',...
                    'Up              -> Increase contour level',...
                    'Down            -> Decrease contour level',...
                    'Left click      -> Select closest contour to keep',...
                    'Right click     -> Delete closest contour',...
                    'C               -> Define crop window for contours',...
                    'D               -> Manually draw points as a contour, left click to draw, right to finish',...
                    'V               -> Add selected contours from the previous slice for use in current slice',...
                    'O               -> Reorder contour points. Warning this is a slow process work in ordered e.g. clockwise fashion to avoid ordering',...
                    'R               -> Remove all selected/drawn contours',...
                    'S               -> Create split contour, i.e. keep current contour segment but add aditional for this slice',...
                    'Space bar       -> Done with current slice, proceed to next',...
                    'H               -> Display help',... 
                    'M               -> Alter colorbar limits',...
                    'L               -> Jump to contour level',...
                    'J               -> Change contour increment size (increment when up key is pressed)',...
                    };
                helpButton = questdlg(msgText,'Help information','OK','OK');
            otherwise
                titleStringAdd=[' Invalid input "',num2str(b),'", press H for help'];                
        end
    end
    
    if ~isempty(Vcs{qSlice}) %Save if not empty
        %Recovery operations
        recoveryStruct.Vcs=Vcs{qSlice};
        recoveryStruct.splitInd=splitInd;
        recoveryStruct.V=V;
        recoveryStruct.C=C;
        recoveryStruct.ic=ic;
        recoveryStruct.ic_old=ic_old;
        recoveryStruct.minC=minC;
        recoveryStruct.maxC=maxC;
        
        save(recoveryFileName,'recoveryStruct');
    end
end

close(hf1);

%% SAVING DATA
if ~isempty(saveName)
    save(saveName,'Vcs');
end

%% DELETING RECOVERY FILES
for qSlice=sliceRange
    recoveryFileName=['recoverFile_uiContourSegment_',num2str(qSlice)];
    if exist(recoveryFileName,'file')==2
        delete(recoveryFileName);
    end
end

end

function setDefaultPointer
   set(gcf,'Pointer','watch'); 
end


