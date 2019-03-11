function [varargout]=hex8_hex20(varargin)

% function [E_hex20,V_hex20,V_hex20_cell,Fb_hex20,S]=hex8_hex20(E_hex8,V_hex8,V_hex8_cell,Fb_hex8)
% ------------------------------------------------------------------------
% This function converts 8 node (e.g. linear) hexahedral elements into 20
% node (e.g. quadratic) hexahedral elements compatible with FEBio. 
%
%
% varargout{1}=E_hex20;
% varargout{2}=V_hex20;
% varargout{3}=V_hex20_cell;
% varargout{4}=Fb_hex20;
% varargout{5}=S;
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/10/19 Created based on tet4_tet10
%------------------------------------------------------------------------
%%

%% Parse input
switch nargin
    case 2
        E_hex8=varargin{1};
        V_hex8=varargin{2};
        V_hex8_cell={};    
        Fb_hex8={};
    case 3
        E_hex8=varargin{1};
        V_hex8=varargin{2};
        V_hex8_cell=varargin{3};
        Fb_hex8={};
    case 4
        E_hex8=varargin{1};
        V_hex8=varargin{2};
        V_hex8_cell=varargin{3};
        Fb_hex8=varargin{4};
    otherwise
        error('wrong number of input arguments');
end

%% Create element set

%Creating "egdge indices" matrices
matVirtSize=size(V_hex8,1)*ones(1,2);

%Edges matrix
indMatrix_1=[E_hex8(:,1) E_hex8(:,2) E_hex8(:,3) E_hex8(:,4)  E_hex8(:,5) E_hex8(:,6) E_hex8(:,7) E_hex8(:,8)  E_hex8(:,5) E_hex8(:,6) E_hex8(:,7) E_hex8(:,8)];
indMatrix_2=[E_hex8(:,2) E_hex8(:,3) E_hex8(:,4) E_hex8(:,1)  E_hex8(:,6) E_hex8(:,7) E_hex8(:,8) E_hex8(:,5)  E_hex8(:,1) E_hex8(:,2) E_hex8(:,3) E_hex8(:,4)];

%3D edges matrix, first layer, first points, second layer, second points
edgeMatrix=indMatrix_1;
edgeMatrix(:,:,2)=indMatrix_2;
edgeMatrix=sort(edgeMatrix,3);

%Convert to virtual indices
edgeIndexMatrix=sub2ind(matVirtSize,edgeMatrix(:,:,1),edgeMatrix(:,:,2));

%Compose unique set
[edgeIndexUni,ind1,ind2]=unique(edgeIndexMatrix(:));

indhex20_new=reshape(ind2,size(edgeIndexMatrix));

%Create unique point index matrices
[indMatrix_1_uni,indMatrix_2_uni]=ind2sub(matVirtSize,edgeIndexUni);

E_hex20=[E_hex8 indhex20_new+size(V_hex8,1)];

%% Create and add new coordinates
% Calculate coordinates for unique new points
V_new=0.5*(V_hex8(indMatrix_1_uni,:)+V_hex8(indMatrix_2_uni,:));%
V_hex20=[V_hex8;V_new];

%% Process cell data 

if ~isempty(V_hex8_cell)
    V_hex20_cell=V_hex8_cell;
    for qc=1:1:numel(V_hex20_cell)
        XX_hex8=V_hex8_cell{qc};        
        XX_new=0.5*(XX_hex8(indMatrix_1_uni,:)+XX_hex8(indMatrix_2_uni,:));
        XX_hex20=[XX_hex8;XX_new];
        V_hex20_cell{qc}=XX_hex20;
    end
else
    V_hex20_cell={};
end


%% Compose output
varargout{1}=E_hex20;
varargout{2}=V_hex20;
varargout{3}=V_hex20_cell;

if nargout>3
    
    % Process boundary faces
    if ~isempty(Fb_hex8)
        
        S = sparse(indMatrix_1_uni,indMatrix_2_uni,ind2(ind1),matVirtSize(1),matVirtSize(2),numel(ind1));
        
        %Edges matrix
        indMatrix_1=[Fb_hex8(:,1) Fb_hex8(:,2) Fb_hex8(:,3) Fb_hex8(:,4)];
        indMatrix_2=[Fb_hex8(:,2) Fb_hex8(:,3) Fb_hex8(:,4) Fb_hex8(:,1)];
        
        %3D edges matrix, first layer, first points, second layer, second points
        edgeMatrix=indMatrix_1;
        edgeMatrix(:,:,2)=indMatrix_2;
        edgeMatrix=sort(edgeMatrix,3);
        
        edgeIndexMatrix=sub2ind(matVirtSize,edgeMatrix(:,:,1),edgeMatrix(:,:,2));
        
        indhex20_new=full(S(edgeIndexMatrix));
        
        Fb_hex20_quad8=[Fb_hex8 indhex20_new+size(V_hex8,1)];
        Fb_hex20=Fb_hex20_quad8(:,[1 5 2 6 3 7 4 8]);
    else
        Fb_hex20_quad8=[];
        Fb_hex20=[];
    end
    varargout{4}=Fb_hex20;
    varargout{5}=Fb_hex20_quad8;
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
