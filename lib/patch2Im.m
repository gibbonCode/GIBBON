function [varargout]=patch2Im(varargin)

% function [M,G,bwLabels]=patch2Im(F,V,C,voxelSize,imOrigin,imSiz)
% -----------------------------------------------------------------------
% This function converts the input surface, specified by the
% faces F and the vertices V into an image. The image size, origin and
% voxel size can be specified. If not specified the image voxel size is set
% to the mean edge length, the size is set to fit the object (with an extra
% voxel spacing all around) and the origin is choosen in the corner of the
% tightly fit image.
% If required the surface is resampled to subvoxel resolution using the
% |subtri| function. Then surface vertices are simply mapped to an image
% coordinate system and since they are densely sampled with respect to the
% voxel size form an enclosing boundary of voxels. The exterior, boundary
% and interior are then formulated using the bwlabeln function. The
% exterior, boundary and intertior voxels are labelled using 0's, 1's and
% 2's in the output image M. The second output G is a structure containing
% the fields G.voxelSize and G.origin which form the geometric information
% for the image.
%
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 22/12/2015 Based on triSurf2Im, expanded for patches in general
%------------------------------------------------------------------------

%%

if nargin<2 || nargin>6
    error('Wrong number of input arguments!');
end

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=[];
        voxelSize=[];
        imOrigin=[];
        imSiz=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=[];
        imOrigin=[];
        imSiz=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=[];
        imSiz=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=varargin{5};
        imSiz=[];
    case 6
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=varargin{5};
        imSiz=varargin{6};
    otherwise
        error('Wrong number of input arguments!');
end

%Removing unused points
[F,V]=removeNotIndexed(F,V);

if isempty(C)
    C=ones(size(F,1),1);
end

if any(C==0)
    error('Zero entries found in C, this label is reserved for boundary voxels, consider using integer labels larger than 0')
end

%% Convert to triangulated surface

nNode=size(F,2);
if nNode==4 %Quadrilateral faces
    [F,V,C]=quad2tri(F,V,'b',C);
elseif nNode>4
    [F,V,C]=patch2tri(F,V,C);
end

%Get colors
C_uni=sort(unique(C(:))); %Unique color set
numC=numel(C_uni); %Number of colors

%Sort order of processing
[~,indSort]=sort(abs(C_uni));
C_uni_sort=C_uni(indSort);

%Compute image and labels for all boundaries
[MT,G,bwLabels]=triSurf2Im(F,V,voxelSize,imOrigin,imSiz);
M=MT;
M(MT==0)=NaN;
M(MT==1)=0;
M(MT==2)=C_uni_sort(1);

if numC>1
    %Store geometry and size information
    imSiz=size(M);
    voxelSize=G.voxelSize;
    imOrigin=G.origin;
    
    for q=2:1:numC;
        [Fs,Vs]=removeNotIndexed(F(C==C_uni_sort(q),:),V);
        [M_s]=triSurf2Im(Fs,Vs,voxelSize,imOrigin,imSiz);
        M(M_s==1)=0; %Boundary voxels
        M(M_s==2)=C_uni_sort(q); %Interior region
    end
end
%%
varargout{1}=M;
varargout{2}=G;
varargout{3}=bwLabels;

 
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
