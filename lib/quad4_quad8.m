function [varargout]=quad4_quad8(varargin)

% function [QUAD8,V8,VX8C]=quad4_quad8(QUAD4,V4,VXC)
%
% This function converts 4 node (e.g. linear) quadrilaterial elements into
% 8 node (e.g. quadratic) quadrilateral elements compatible with FEBio.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/10/19 Created based on tri3_tri6
%----------------------------------------------------------------------
%%

switch nargin
    case 2
        QUAD4=varargin{1};
        V4=varargin{2};
        VXC={};
    case 3
        QUAD4=varargin{1};
        V4=varargin{2};
        VXC=varargin{3};
end

%%

E=[QUAD4(:,[1 2]); QUAD4(:,[2 3]); QUAD4(:,[3 4]); QUAD4(:,[4 1])]; %Edges matrix
Es=sort(E,2); %Sorted edges matrix
[~,ind1,~]=unique(Es,'rows'); %Indices for unique edges
E=E(ind1,:); %The unique esges

numPoints = size(V4,1);
numEdges = size(E,1);

% Get indices of the three edges associated with each face
A = sparse(E(:,1),E(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
A = max(A,A'); %Copy symmetric

%Indices for A matrix
indA_12=QUAD4(:,1)+(QUAD4(:,2)-1)*numPoints;
indA_23=QUAD4(:,2)+(QUAD4(:,3)-1)*numPoints;
indA_34=QUAD4(:,3)+(QUAD4(:,4)-1)*numPoints;
indA_41=QUAD4(:,4)+(QUAD4(:,1)-1)*numPoints;

%Get indices for vertex array
indV_12=full(A(indA_12));
indV_23=full(A(indA_23));
indV_34=full(A(indA_34));
indV_41=full(A(indA_41));

%Create faces array
QUAD8=[QUAD4(:,1) QUAD4(:,2) QUAD4(:,3) QUAD4(:,4) indV_12 indV_23 indV_34 indV_41];

%Create vertex array
Vn=0.5*(V4(E(:,1),:)+V4(E(:,2),:)); %new mid-edge points
V8 = [V4; Vn]; %Join point sets

%%
varargout{1}=QUAD8;
varargout{2}=V8;

if nargout==3
    %Derive VX8C
    if ~isempty(VXC)
        VX8C=VXC;
        for q=1:1:numel(VXC)
            VX=VXC{q};
            VX_1_4=VX;
            VX_5_8=0.5*(VX(E(:,1),:)+VX(E(:,2),:)); %new mid-edge data
            VX8=[VX_1_4; VX_5_8];
            VX8C{q}=VX8;
        end
    else
        VX8C={};
    end
    varargout{3}=VX8C;
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
