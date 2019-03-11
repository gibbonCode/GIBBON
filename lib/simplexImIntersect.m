function [regionLabel]=simplexImIntersect(varargin)

% [logicOut]=simplexImIntersect(F1,V1,C1,V2,voxelSize)
% ------------------------------------------------------------------------
% This function uses the input simplex (or simplices) defined by F1
% (faces), V1 (vertices), and C1 (boundary labels) to label the vertices V2
% according to what simplex regions they are contained in. The labels are
% stored in the output regionLabel. A vertex inside a simplex region is
% labelled with the boundary label of the simplex region. Vertices outside
% of all simplices are labelled with NaN. Points are labelled based on the
% patch2Im function, i.e. the following steps are used: 
% 1) The simplex is converted to an image where voxels intensities denote
% simplex region labels. (see patch2Im).
% 2) The vertices in V2 are converted to image coordinates in this image
% 3) The image coordinates are converted to image voxel indices
% 4) Voxel indices are used to retrieve the labels from the simplex image. 
%
% The image constructed uses the optional input voxelSize (the default if
% not provided is half of the mean edge size of the input simplex). The
% completeness/accuracy of the labelling depends on the voxel size. Some
% points are labelled as 0 which means they are found inside boundary
% voxels and cannot, based on the current voxel size, be assigned as
% outside or inside a particular region. All vertices labelled as NaN or a
% value >0 are labelled correctly. However some of the vertices labelled as
% 0 may actually be inside a simplex region or outside all simplex regions,
% i.e. 0 denotes that their status is unknown given the voxel size used. 
%
% See also |patch2im|
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/02/07:  Created
% 2018/02/08:  Expanded to allow for labeling and use with multi-region
% simplices. 
%------------------------------------------------------------------------

%% Parse input 

switch nargin    
    case 4
        F1=varargin{1};
        V1=varargin{2};
        C1=varargin{3};
        V2=varargin{4};
        voxelSize=[];
    case 5
        F1=varargin{1};
        V1=varargin{2};
        C1=varargin{3};
        V2=varargin{4};
        voxelSize=varargin{5};
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
[M]=patch2Im(F1,V1,C1,voxelSize,imOrigin,imSiz);
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
regionLabel=nan(size(V2,1),1); %Initialize as false
regionLabel(logicValid)=M(indValid);

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
