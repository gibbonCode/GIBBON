function [varargout]=patch2Im(varargin)

% function [M,G,bwLabels]=patch2Im(F,V,C,voxelSize,imOrigin,imSiz,boundaryType)
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
% 2015/12/22 Based on triSurf2Im, expanded for patches in general
% 2020/11/12 Adding handling of boundary types
%------------------------------------------------------------------------

%%

if nargin<2 || nargin>7
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
        boundaryType=0;
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=[];
        imOrigin=[];
        imSiz=[];
        boundaryType=0;
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=[];
        imSiz=[];
        boundaryType=0;
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=varargin{5};
        imSiz=[];
        boundaryType=0;
    case 6
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=varargin{5};
        imSiz=varargin{6};
        boundaryType=0;
    case 7 
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        voxelSize=varargin{4};
        imOrigin=varargin{5};
        imSiz=varargin{6};
        boundaryType=varargin{7};
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

%%
%Compute image and labels for all boundaries
[MT,G,bwLabels]=triSurf2Im(F,V,voxelSize,imOrigin,imSiz);

%Store geometry and size information
imSiz=size(MT);
voxelSize=G.voxelSize;
imOrigin=G.origin;
M=NaN(imSiz); %Initialize as all NaN

if numC>1
    for q=1:1:numC
        [Fs,Vs]=removeNotIndexed(F(C==C_uni_sort(q),:),V);
        [M_s]=triSurf2Im(Fs,Vs,voxelSize,imOrigin,imSiz);
        switch boundaryType
            case -1 %Exclusive, regard boundary as out
                M(M_s==2)=C_uni_sort(q); %Interior region
            case 0
                M(M_s==1)=0; %Boundary voxels
                M(M_s==2)=C_uni_sort(q); %Interior region
            case 1 %Inclusive, regard boundary as in
                M(M_s>0)=C_uni_sort(q); %Interior region
        end
    end
else
    switch boundaryType
        case -1 %Exclusive, regard boundary as out
            M(MT==2)=C_uni_sort; %Interior region
        case 0
            M(MT==1)=0; %Boundary voxels
            M(MT==2)=C_uni_sort; %Interior region
        case 1 %Inclusive, regard boundary as in
            M(MT>0)=C_uni_sort; %Interior region
    end
end

%% Gather output
varargout{1}=M;
varargout{2}=G;
varargout{3}=bwLabels;

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
