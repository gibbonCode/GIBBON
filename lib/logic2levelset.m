function [D]=logic2levelset(varargin)

% function [D]=logic2levelset(logicInside,voxelSize,voxelSizeResample)
% ------------------------------------------------------------------------
% This function converts the logic image logicInside to a level set image.
% The level set image intensities define the distance of the voxels to the
% boundaries of the logic image. 
% 
% 
% Change log: 
% 
% ------------------------------------------------------------------------
%%

switch nargin
    case 1
        logicInside=varargin{1};
        voxelSize=[];
        voxelSizeResample=[];
    case 2
        logicInside=varargin{1};
        voxelSize=varargin{2};
        voxelSizeResample=[];
    case 3
        logicInside=varargin{1};
        voxelSize=varargin{2};
        voxelSizeResample=varargin{3};
end

if isempty(voxelSize)
    voxelSize=ones(1,3);
end

if isempty(voxelSizeResample)
    voxelSizeResample=mean(voxelSize);
end

if numel(voxelSizeResample)==1
    voxelSizeResample=voxelSizeResample*ones(1,3);
end

%%

siz=size(logicInside);
%Resample input image isotropically (isotropic voxels)
if ~max(voxelSize-voxelSizeResample)< max(eps(voxelSize))
    resampleOn=1;    
    [logicInside]=imageResample(logicInside,voxelSize,voxelSizeResample);
    logicInside=logicInside>0; %Forces logic and fixes potential NaN's
    scaleFactor=voxelSizeResample(1);
else
    resampleOn=0;
    scaleFactor=voxelSize(1);
end
logicInside=logicInside>0; %Force to be a logic (in case resampling altered it)

%Remove interior from logic
try    
    logicOn = bwmorph3(logicInside,'remove'); %New in R2018a
catch
    logicOn=logicRemoveInterior(logicInside); %GIBBON alterative
end

%Do distance transform on isotropic image
D = double(bwdist(logicOn,'euclidean')); %Compute distance
D(logicInside)=-D(logicInside); %Negate distance inside
D=D.*scaleFactor; %Scale by voxel size

if resampleOn==1    
    sizNew=size(D); %Size of distance data in resampled state
    [J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3)); %Image coordinate grid of original image
    [X,Y,Z]=im2cart(I,J,K,voxelSize); %Spatial coordinate of original image grid
    [I,J,K]=cart2im(X,Y,Z,voxelSizeResample); %Image coordinates of original grid in resampled image
    I=round(I); J=round(J); K=round(K); %Rounding 
    IND=reshape(sub2indn(sizNew,[I(:) J(:) K(:)],1),size(I)); %Convert to linear indices
    D=D(IND); %Override D to be data sampled at original image points
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
