function [varargout]=subTriDual(varargin)

% function [FT,VT,C_type,indIni,C_new]=subTriDual(F,V,logicFaces,C)
% ------------------------------------------------------------------------
% This function can globally or locally refine the input triangulation
% based on so-called "dual refinement". The process of refinement can be
% visualized as first taking the dual, and to triangulate the dual by
% incorporating the original points as well. 
%
% The input consists of:
% F:          the faces
% V:          the vertices
% logicFaces: A logic for the faces requiring refinement
% C:          color data on either the faces or the vertices
% 
% The ouput can consist of: 
% FT:     Faces
% VT:     Vertices
% C_type: Color/label for triangle type, i.e. original (1), refined (2), or
%         boundary (3) 
% indIni: Indices for original points
% C_new:  New color data for faces or vertices
%
% Change log: 
% 2010/01/01 Added for GIBBON
% 2019/06/12 Fixed bug when logic contains "sharp teeth"
% 2019/06/12 Rewrote whole function in a different way using most recent
% tools and capabilities. Code is now much more efficient, handles face or
% vertex color data, fixes several bugs in old code. 
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        logicFaces=true(size(F,1),1);
        C=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        logicFaces=varargin{3};
        C=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        logicFaces=varargin{3};
        C=varargin{4};
end

if isempty(logicFaces)
    logicFaces=true(size(F,1),1);
end

%%

% Get connectivity data
conStruct=patchConnectivity(F(logicFaces,:),V,{'ef','ev'}); %Connectivity structure
EF=conStruct.edge.face; %Edge-face connectivity
E=conStruct.edge.vertex; %Edge-vertex connectivity

% Check boundary regions
logicBoundary=any(EF==0,2);
indF=EF(logicBoundary,1);

% Create new coordinate array
VF=patchCentre(F(logicFaces,:),V); %New face centre coordinates
VT=[V;VF]; %The new vertex array

% Compone new face matrix

e=E(~logicBoundary,:); %Edge-vertex indices not part of boundary
ef=EF(~logicBoundary,:); %Edge-face indices not part of boundary

FT=[F(~logicFaces,:);... %The unaltered faces
    e(:,1) ef(:,1)+size(V,1) ef(:,2)+size(V,1);... %New refined face 1 
    e(:,2) ef(:,2)+size(V,1) ef(:,1)+size(V,1);... %New refined face 2
    E(logicBoundary,:) indF+size(V,1);... %New boundary face
    ];

% Get color/label data for original (1), refined (2), or boundary (3)
C_type=[ones(nnz(~logicFaces),1); 2*ones(2*nnz(~logicBoundary),1);3*ones(nnz(logicBoundary),1); ];

%Derive color data for refined set
if ~isempty(C) && nargout>4
    if size(C,1)==size(F,1) %Face color data
        indFaces=find(logicFaces);
        C_new=[C(~logicFaces,:); repmat((C(indFaces(ef(:,1)),:)+C(indFaces(ef(:,2)),:))/2,2,1); C(indFaces(indF),:); ];
    elseif size(C,1)==size(V,1) %Vertex color data
        CC=vertexToFaceMeasure(F(logicFaces,:),C);
        C_new=[C;CC];
    else 
        error('Color data should be nxq in size whereby n is the number of faces or the number of vertices');
    end
else
    C_new=[];
end

%% Collect output
varargout{1}=FT; %Faces
varargout{2}=VT; %Vertices
varargout{3}=C_type; %Color/label for triangle type
varargout{4}=1:1:size(V,1); %Indices for original points
varargout{5}=C_new; %New color data for faces or vertices

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
