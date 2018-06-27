function [varargout]=tet4_tet10(varargin)

% function [E_tet10,V_tet10,V_tet10_cell,Fb_tet10,S]=tet4_tet10(E_tet4,V_tet4,V_tet4_cell,Fb_tet4)
% ------------------------------------------------------------------------
% This function converts 4 node (e.g. linear) tetrahedral elements into 10
% node (e.g. quadratic) tetrahedral elements compatible with FEBio. 
%
%
% varargout{1}=E_tet10;
% varargout{2}=V_tet10;
% varargout{3}=V_tet10_cell;
% varargout{4}=Fb_tet10;
% varargout{5}=S;
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
% 2015/09/22 Renamed variables
% 2015/09/22 Updated "bookkeeping" of indices of added points to avoid
% unique operation
% 2015/09/22 Added "bookkeeping" of faces e.g. boundaries
% 2015/09/22 Added sparse index array output option
%------------------------------------------------------------------------
%%

%% Parse input
switch nargin
    case 2
        E_tet4=varargin{1};
        V_tet4=varargin{2};
        V_tet4_cell={};    
        Fb_tet4={};
    case 3
        E_tet4=varargin{1};
        V_tet4=varargin{2};
        V_tet4_cell=varargin{3};
        Fb_tet4={};
    case 4
        E_tet4=varargin{1};
        V_tet4=varargin{2};
        V_tet4_cell=varargin{3};
        Fb_tet4=varargin{4};
    otherwise
        error('wrong number of input arguments');
end

%%

%Creating "egdge indices" matrices
matVirtSize=size(V_tet4,1)*ones(1,2);

%Edges matrix
indMatrix_1=[E_tet4(:,1) E_tet4(:,2) E_tet4(:,3) E_tet4(:,1) E_tet4(:,2) E_tet4(:,3)];
indMatrix_2=[E_tet4(:,2) E_tet4(:,3) E_tet4(:,1) E_tet4(:,4) E_tet4(:,4) E_tet4(:,4)];

%3D edges matrix, first layer, first points, second layer, second points
edgeMatrix=indMatrix_1;
edgeMatrix(:,:,2)=indMatrix_2;
edgeMatrix=sort(edgeMatrix,3);

%Convert to virtual indices
edgeIndexMatrix=sub2ind(matVirtSize,edgeMatrix(:,:,1),edgeMatrix(:,:,2));

%Compose unique set
[edgeIndexUni,ind1,ind2]=unique(edgeIndexMatrix(:));

indTet10_new=reshape(ind2,size(edgeIndexMatrix));

%Create unique point index matrices
[indMatrix_1_uni,indMatrix_2_uni]=ind2sub(matVirtSize,edgeIndexUni);

S = sparse(indMatrix_1_uni,indMatrix_2_uni,ind2(ind1),matVirtSize(1),matVirtSize(2),numel(ind1));

E_tet10=[E_tet4 indTet10_new+size(V_tet4,1)];


%%
%Calculate coordinates for unique new points
V_tet10_cell=V_tet4_cell;

V_new=zeros(numel(edgeIndexUni),size(V_tet4,2));
for q=1:1:size(V_tet4,2)
    X=V_tet4(:,q);
    Xq=0.5*(X(indMatrix_1_uni)+X(indMatrix_2_uni));%
    V_new(:,q)=Xq;
end
V_tet10=[V_tet4;V_new];


if ~isempty(V_tet4_cell)
    for qc=1:1:numel(V_tet10_cell)
        XX_tet4=V_tet4_cell{qc};
        XX_new=zeros(numel(edgeIndexUni),size(XX_tet4,2));
        for q=1:1:size(XX_tet4,2)
            X=XX_tet4(:,q);
            Xq=0.5*(X(indMatrix_1_uni)+X(indMatrix_2_uni));%
            XX_new(:,q)=Xq;
        end
        XX_tet10=[XX_tet4;XX_new];
        V_tet10_cell{qc}=XX_tet10;
    end
end

%%

if ~isempty(Fb_tet4)
    
    %Edges matrix
    indMatrix_1=[Fb_tet4(:,1) Fb_tet4(:,2) Fb_tet4(:,3)];
    indMatrix_2=[Fb_tet4(:,2) Fb_tet4(:,3) Fb_tet4(:,1)];
    
    %3D edges matrix, first layer, first points, second layer, second points
    edgeMatrix=indMatrix_1;
    edgeMatrix(:,:,2)=indMatrix_2;
    edgeMatrix=sort(edgeMatrix,3);
    
    edgeIndexMatrix=sub2ind(matVirtSize,edgeMatrix(:,:,1),edgeMatrix(:,:,2));
    
    indTet10_new=full(S(edgeIndexMatrix));
    
    Fb_tet10=[Fb_tet4 indTet10_new+size(V_tet4,1)];
    Fb_tet10=Fb_tet10(:,[1 4 2 5 3 6]);
else
    Fb_tet10=[];
end

%% Compose output
varargout{1}=E_tet10;
varargout{2}=V_tet10;
varargout{3}=V_tet10_cell;
varargout{4}=Fb_tet10;
varargout{5}=S;   
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
