function [logicOut]=simplexImIntersect(varargin)

% [logicOut]=simplexImIntersect(F1,V1,V2,voxelSize)
% ------------------------------------------------------------------------
% This function computes a logic which is true for vertices in V2 which are
% deemed outside of the simplex defined by the faces array F1 and the
% vertices array V1. Points are deemed outside based on the patch2Im
% function, i.e. the following steps are used: 
% 1) The simplex is converted to an image where voxels are 1 if they are in
% or on the simplex and 0 if the are outside. 
% 2) The vertices in V2 are converted to image coordinates in this image
% 3) The image coordinates are converted to image voxel indices
% 4) Voxel indices which point at voxels which are 0 in the simplex image
% are outside. 
%
% The image constructed uses the optional input voxelSize (default if not
% provided is half of the mean edge size of the input simplex). Although
% all of the vertices where logicOut is 1 are out not all that are 0 in
% logicOut are truely in. The level of accuracy of the logic depends on the
% voxel size. If set too high then too many points are deemed outside. If
% set too small then computational time is very large. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/02/07:  Created
%
%------------------------------------------------------------------------

%% Parse input 

switch nargin    
    case 3
        F1=varargin{1};
        V1=varargin{2};
        V2=varargin{3};
        voxelSize=[];
    case 4
        F1=varargin{1};
        V1=varargin{2};
        V2=varargin{3};
        voxelSize=varargin{4};
end

if isempty(voxelSize)
    [D1]=patchEdgeLengths(F1,V1);
    voxelSize=mean(D1(:))/2;
end

%% Compute simplex image
% Using |patch2Im| function to convert patch data to image data

% Settings 
imOrigin=min(V1,[],1)-voxelSize;
imMax=max(V1,[],1)+voxelSize;
imSiz=round((imMax-imOrigin)/voxelSize);
imSiz=imSiz([2 1 3]); %Image size (x, y corresponds to j,i in image coordinates, hence the permutation)

% Compute image
[M]=patch2Im(F1,V1,[],voxelSize,imOrigin,imSiz);
logicImage=~isnan(M); %Logic image for voxels in or on the simplex

%% Find vertices of simplex 2 outside of simplex 1

indImage=find(logicImage); %Linear indices for voxels in or on the simplex

% Convert to image coordinates
V2_im=V2-imOrigin(ones(size(V2,1),1),:);
[V2_im(:,1),V2_im(:,2),V2_im(:,3)]=cart2im(V2_im(:,1),V2_im(:,2),V2_im(:,3),voxelSize*ones(1,3));
V2_im=round(V2_im);

logicValid=all(V2_im>1,2) & V2_im(:,1)<=size(M,1) & V2_im(:,2)<=size(M,2) & V2_im(:,3)<=size(M,3); %Valid voxels only have indices larger than 1

indValid = sub2ind(size(M),V2_im(logicValid,1),V2_im(logicValid,2),V2_im(logicValid,3));

%Calculate the outside logic
logicOut=false(size(V2,1),1); %Initialize as false
logicOut(logicValid)=~ismember(indValid,indImage); %Check if indices are a member of the in or on set
logicOut=logicOut | ~logicValid; %Add values that were already deemed outside

