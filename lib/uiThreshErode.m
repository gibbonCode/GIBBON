function [L_BG]=uiThreshErode(varargin)

% function [L_BG]=uiThreshErode(M,thresholdInitial,blurKernelSize,groupCropOption,inConvexHull)


%% Parse input

thresholdInitial=[];
blurKernelSize=[];
groupCropOpt=[];
inConvexHullOpt=[];
switch nargin
    case 1
        M=varargin{1};
    case 2
        M=varargin{1};
        thresholdInitial=varargin{2};
    case 3
        M=varargin{1};
        thresholdInitial=varargin{2};
        blurKernelSize=varargin{3};
    case 4
        M=varargin{1};
        thresholdInitial=varargin{2};
        blurKernelSize=varargin{3};
        groupCropOpt=varargin{4};
    case 5
        M=varargin{1};
        thresholdInitial=varargin{2};
        blurKernelSize=varargin{3};
        groupCropOpt=varargin{4};
        inConvexHullOpt=varargin{5};
end

%Set defaults
if isempty(thresholdInitial)
     thresholdInitial=0.1;
end

%Quick fix for low thresholds 
if thresholdInitial<0.01
    thresholdInitial=0.01;
end

if isempty(blurKernelSize)
     blurKernelSize=ndims(M);
end

if isempty(groupCropOpt)
     groupCropOpt=0;
end

if isempty(inConvexHullOpt)
     inConvexHullOpt=0;
end
%% PLOT SETTINGS

fontSize=10;
cMap=gray(250);
falpha=1;

logicVoxels=false(size(M));
logicVoxels(round(size(M,1)/2),:,:)=1;
logicVoxels(:,round(size(M,2)/2),:)=1;
logicVoxels(:,:,round(size(M,3)/2))=1;

%%
M_original=M;
M=double(M);

%% Normalizing the image data
M=dataNorm(M,1);

%% APPLY PRE-BLURRING
if blurKernelSize>0
    %Defining erosion/dilation kernel
    k=blurKernelSize;
    p=k-round(k./2);
    
    hb=gauss_kernel(k,ndims(M),1.5,'width'); %Gaussian kernel
       
    M_rep=zeros(size(M)+(2.*p));
    M_rep(p+(1:size(M,1)),p+(1:size(M,2)),p+(1:size(M,3)))=M;
    M = convn(M_rep,hb,'valid');
end

%% Defining erosion/dilation kernel
k=3;
hg=ones(k,k,k); 
p=k-round(k./2);

%% Start loop
done=0;
while done==0
    %Get initial threshold level    
    T_threshold=thresholdInitial; %default start value
    done=0;
    qc=1;
    runMode=1;
    while done==0
        
        switch runMode
            case 1 %Determining threshold
                L_BG=M>T_threshold;
            case 2 %Setting erosion/dilation
                %"Blurring current mask"
                L_BG_rep=zeros(size(M)+(2.*p));
                L_BG_rep(p+(1:size(M,1)),p+(1:size(M,2)),p+(1:size(M,3)))=L_BG;
                L_BG_blur = convn(double(L_BG_rep),hg,'valid');
        end
        
        if qc==1 %Plot figure on first iteration
            
            %Cropped image            
            [Fs,Vs,Cs]=ind2patch(logicVoxels&L_BG,M_original,'v');
            
            hf1=cFigure;
            title(['Threshold is ',num2str(T_threshold),'*mean, press up to increase or down to decrease (by 10%), press space to keep and continue'],'FontSize',fontSize);
            hold on; xlabel('X-J','FontSize',fontSize);ylabel('Y-I','FontSize',fontSize);zlabel('Z-K','FontSize',fontSize);
            hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat');
            set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); colorbar; 
            set(gca,'FontSize',fontSize);
            drawnow;
            
        else %redefine patch data
            delete(hs); %remove patch data from figure
            
            %redefine patch data for the cropped image
            [Fs,Vs,Cs]=ind2patch(logicVoxels&L_BG,M_original,'v');
            
            %patch again
            if runMode==1
                title(['Threshold is ',num2str(T_threshold),'*mean, press up to increase or down to decrease (by 10%), press space to keep and continue'],'FontSize',fontSize);
            else
                title('Current cropping press up/down to dilate/erode, WARNING EROSION MAY REMOVE ENTIRE SLICES','FontSize',fontSize);
            end
            
            hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat');
            set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); colorbar; 
            set(gca,'FontSize',fontSize);
            drawnow;
        end
        
        [~,~,b]=qginput(1);
        switch b
            case 30 %  up arrow key => increase threshold by 10% OR dilate
                if runMode==1
                    T_threshold=T_threshold.*1.1;
                else
                    L_BG=(L_BG_blur>0); %blurred boundary elements that increased due to blur are now set to 1 as well so added to mask
                end
            case 31 %down arrow key => decrease threshold by 10% OR erode
                if runMode==1
                    T_threshold=T_threshold.*0.9;
                else
                    L_BG=(L_BG_blur==sum(hg(:))); %blurred boundary elements that decreased due to blur are now set to 0 as well so removed from mask
                end
            case 32 %spacebar key => confirms current threshold result, exits while loop
                if runMode==1
                    runMode=2;
                else
                    done=1;
                end
        end
        qc=qc+1;
    end
    
    %% GROUPING BASED CROPPING
    if groupCropOpt==1
        
        GROUP_STRUCT = bwconncomp(L_BG,6);
        IMAGE_OBJECTS=GROUP_STRUCT.PixelIdxList;
        
        %Filtering out groups that are too small
        group_sizes = cell2mat(cellfun(@(x) max(size(x)), IMAGE_OBJECTS,'UniformOutput',0)'); %get group sizes
        [~,ind_max]=max(group_sizes); %Index of largest group
        IND_OBJECT=IMAGE_OBJECTS{ind_max}; %indices of voxels for largest group
        
        %Redefining background logic using group
        L_BG=false(size(L_BG)); 
        L_BG(IND_OBJECT)=1;   
        
        %Plotting the result
        delete(hs); %remove patch data from figure
        
        %redefine patch data for the cropped image
        [Fs,Vs,Cs]=ind2patch(logicVoxels&L_BG,M_original,'v');
                
        title('Grouping result','FontSize',fontSize);
               
        hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat');
        set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); 
        set(gca,'FontSize',fontSize);
        drawnow;
    end

    %% Convex hull based filling
    if inConvexHullOpt==1
        [I,J]=ndgrid(1:1:size(M,1),1:1:size(M,2));
        for q=1:1:size(M,3);           
            l_bg=L_BG(:,:,q);            
            if nnz(l_bg)>0
                if nnz(l_bg)>3
                    try
                        P=[I(l_bg) J(l_bg)];
                        K = convhull(P(:,1),P(:,2));
                        indNot=find(~l_bg);
                        [Lp,~] = inpolygon(I(indNot),J(indNot),P(K,1),P(K,2));
                        indAdd=indNot(Lp);
                        L=false(size(M,1),size(M,2));
                        L(indAdd)=1;
                        L_BG(:,:,q)=L_BG(:,:,q) | L;
                    end
                end
            else
               L_BG(:,:,q)=0; 
            end
        end
        
        title('Convexhull fill result','FontSize',fontSize);
        delete(hs); %remove patch data from figure
        %redefine patch data for the cropped image
        [Fs,Vs,Cs]=ind2patch(logicVoxels&L_BG,M_original,'v');

        hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',Cs,'FaceColor','flat');
        set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap);
        set(gca,'FontSize',fontSize);
        drawnow;
    end
    
    %%
    Q = questdlg('Would you like to keep this result?','Keep result?','YES','NO','YES');
    switch Q
        case 'YES'
            done=1;
        case 'NO'
            done=0;
    end
    close(hf1)
    
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
