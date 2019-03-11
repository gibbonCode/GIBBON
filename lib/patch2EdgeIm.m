function [varargout]=patch2EdgeIm(varargin)

% function [M,G]=patch2EdgeIm(F,V,voxelSize,imOrigin,imSiz,n)
% -----------------------------------------------------------------------
% This function converts the input patch data specified by the faces F and
% the vertices V into an edge image. 
% The image size, origin, voxel size and the number of edge subdevissions
% can be specified. If not specified the image voxel size is set to 1/10th
% the mean edge length, the size is set to fit the object (with an extra
% voxel spacing all around), the origin is choosen in the corner of the
% tightly fit image and the number of subdevisions is chosen as the maximum
% edge length devided by 1/5th of the voxel size. 
% The output consists of an image M where ones denote edge voxels. The
% optional second output G is a structure containing the fields G.voxelSize
% and G.origin which form the geometric information for the image. 
% 
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/08/2013
%------------------------------------------------------------------------

if nargin<2 || nargin>6
        error('Wrong number of input arguments!');
end

%Get faces and vertices
F=varargin{1};
V=varargin{2};

%Removing unused points
[F,V]=removeNotIndexed(F,V);

%Get edges
E=patchEdges(F,1);

%Checking edge lenghts of surface
[edgeLengths]=patchEdgeLengths(F,V);

maxEdgeLength=max(edgeLengths(:));
meanEdgeLength=mean(edgeLengths(:));

if nargin >2
    voxelSize=varargin{3};
else
    voxelSize=meanEdgeLength/10;
end

if nargin >3
    imOrigin=varargin{4};
else
    %Determine surface set coordinate minima
    minV=min(V,[],1);
    
    %Determine shift so all coordinates are positive
    imOrigin=(minV-voxelSize);
end
    
%Get edge subdevision metric
if nargin >5
    n=varargin{6};
else
    subFactor=5;
    n=maxEdgeLength/(voxelSize/subFactor);
end

%%

% Subdevide edges if required
if (n-1)>eps(1)
    n=ceil(n);
    Vn=zeros(n*size(E,1),size(V,2));
    for q=1:1:size(V,2)
        X=V(:,q);
        XE=X(E);
        Xn=linspacen(XE(:,1),XE(:,2),n);        
        Vn(:,q)=Xn(:);
    end
end

V=Vn;

%Shift points using origin
V=V-imOrigin(ones(size(V,1),1),:);

%Convert to image coordinates
V_IJK=V;
[V_IJK(:,1),V_IJK(:,2),V_IJK(:,3)]=cart2im(V(:,1),V(:,2),V(:,3),voxelSize*ones(1,3));

%Rounding image coordinates to snap to voxel
V_IJK=round(V_IJK);

%Determine image size if not provided
if nargin<5    
    imSiz=max(V_IJK,[],1)+1;
else
    imSiz=varargin{5};
    
    %Remove invalid indices
    V_IJK=V_IJK(V_IJK(:,1)<=imSiz(1) & V_IJK(:,1)>=1,:);
    V_IJK=V_IJK(V_IJK(:,2)<=imSiz(2) & V_IJK(:,2)>=1,:);
    V_IJK=V_IJK(V_IJK(:,3)<=imSiz(3) & V_IJK(:,3)>=1,:);       
end

%Get linear indices of points
indVertices=sub2ind(imSiz,V_IJK(:,1),V_IJK(:,2),V_IJK(:,3));

%Create surface boundary image
M=zeros(imSiz);
M(indVertices)=1;

%Storing image geometry metrics
G.voxelSize=voxelSize;
G.origin=imOrigin;

%Create output
varargout{1}=M;
varargout{2}=G;
    
 
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
