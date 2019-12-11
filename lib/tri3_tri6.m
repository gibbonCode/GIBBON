function [varargout]=tri3_tri6(varargin)

% function [TRI6,V6,VX6C]=tri3_tri6(TRI3,V3,VXC)
%
% This function converts 3 node (e.g. linear) triangular elements into 6
% node (e.g. quadratic) triangular elements compatible with FEBio.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2013/09/12
% 2018/06/26
% 2018/10/19 Altered so merging nodes is no longer required 
%----------------------------------------------------------------------
%%

switch nargin
    case 2
        TRI3=varargin{1};
        V3=varargin{2};
        VXC={};
    case 3
        TRI3=varargin{1};
        V3=varargin{2};
        VXC=varargin{3};
end

%%

E=[TRI3(:,[1 2]); TRI3(:,[2 3]);  TRI3(:,[3 1])]; %Edges matrix
Es=sort(E,2); %Sorted edges matrix
[~,ind1,~]=unique(Es,'rows');
E=E(ind1,:);

numPoints = size(V3,1);
numEdges = size(E,1);

% Get indices of the three edges associated with each face
A = sparse(E(:,1),E(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
A = max(A,A'); %Copy symmetric

%Indices for A matrix
indA_12=TRI3(:,1)+(TRI3(:,2)-1)*numPoints;
indA_23=TRI3(:,2)+(TRI3(:,3)-1)*numPoints;
indA_31=TRI3(:,3)+(TRI3(:,1)-1)*numPoints;

%Get indices for vertex array
indV_12=full(A(indA_12));
indV_23=full(A(indA_23));
indV_31=full(A(indA_31));

%Create faces array
TRI6=[TRI3(:,1) TRI3(:,2) TRI3(:,3) indV_12 indV_23 indV_31];

%Create vertex array
Vn=0.5*(V3(E(:,1),:)+V3(E(:,2),:)); %new mid-edge points
V6 = [V3; Vn]; %Join point sets

%%
varargout{1}=TRI6;
varargout{2}=V6;

if nargout==3
    %Derive VX6C
    if ~isempty(VXC)
        VX6C=VXC;
        for q=1:1:numel(VXC)
            VX=VXC{q};
            VX_1_3=VX;
            VX_4_6=0.5*(VX(E(:,1),:)+VX(E(:,2),:)); %new mid-edge data
            VX6=[VX_1_3; VX_4_6];
            VX6C{q}=VX6;
        end
    else
        VX6C={};
    end
    varargout{3}=VX6C;
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
