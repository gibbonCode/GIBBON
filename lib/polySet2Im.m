function [varargout]=polySet2Im(varargin)

% function [M,G]=polySet2Im(E,V,voxelSize,imOrigin,imSiz,n)
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
E=varargin{1};
V=varargin{2};

%Checking edge lenghts (point spacings) 
[edgeLengths]=patchEdgeLengths(E,V);

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
    siz=max(V_IJK,[],1)+1;
else
    siz=varargin{5};
    
    %Remove invalid indices
    V_IJK=V_IJK(V_IJK(:,1)<=siz(1) & V_IJK(:,1)>=1,:);
    V_IJK=V_IJK(V_IJK(:,2)<=siz(2) & V_IJK(:,2)>=1,:);
    V_IJK=V_IJK(V_IJK(:,3)<=siz(3) & V_IJK(:,3)>=1,:);       
end

%Get linear indices of points
indVertices=sub2ind(siz,V_IJK(:,1),V_IJK(:,2),V_IJK(:,3));

%Create surface boundary image
M=zeros(siz);
M(indVertices)=1;

%Storing image geometry metrics
G.voxelSize=voxelSize;
G.origin=imOrigin;

%Create output
varargout{1}=M;
varargout{2}=G;
    
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
