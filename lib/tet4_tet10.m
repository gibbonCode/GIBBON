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

%% OLD
% 
% %Collect nodes
% V_1_4=V_tet4;
% V_5 =0.5.*(V_tet4(E_tet4(:,1),:)+V_tet4(E_tet4(:,2),:));
% V_6 =0.5.*(V_tet4(E_tet4(:,2),:)+V_tet4(E_tet4(:,3),:));
% V_7 =0.5.*(V_tet4(E_tet4(:,3),:)+V_tet4(E_tet4(:,1),:));
% V_8 =0.5.*(V_tet4(E_tet4(:,1),:)+V_tet4(E_tet4(:,4),:));
% V_9 =0.5.*(V_tet4(E_tet4(:,2),:)+V_tet4(E_tet4(:,4),:));
% V_10=0.5.*(V_tet4(E_tet4(:,3),:)+V_tet4(E_tet4(:,4),:));
% V_tet10=[V_1_4;V_5;V_6;V_7;V_8;V_9;V_10];
% 
% %Define elements
% numTets=size(E_tet4,1);
% num_V_1_4=size(V_1_4,1);
% indTets=(1:numTets)';
% E_tet10=[E_tet4...                            % 1-4
%        indTets+num_V_1_4+(numTets*(1-1))... % 5
%        indTets+num_V_1_4+(numTets*(2-1))... % 6
%        indTets+num_V_1_4+(numTets*(3-1))... % 7
%        indTets+num_V_1_4+(numTets*(4-1))... % 8
%        indTets+num_V_1_4+(numTets*(5-1))... % 9
%        indTets+num_V_1_4+(numTets*(6-1))];  % 10
%    
% %Removing double points
% [~,ind_uni_1,ind_uni_2]=unique(pround(V_tet10,5),'rows');
% V_tet10=V_tet10(ind_uni_1,:);
% E_tet10=ind_uni_2(E_tet10); %Changing indices in faces matrix
% 
% if size(E_tet4,1)==1 %Transpose in this case due to MATLAB behaviour
%     E_tet10=E_tet10';
% end
% 
% %Derive V_tet10_cell
% if ~isempty(V_tet4_cell)    
%     V_tet10_cell=V_tet4_cell;
%     for q=1:1:numel(V_tet4_cell)
%         VX=double(V_tet4_cell{q});
%         VX_1_4=VX;
%         VX_5 =0.5.*(VX(E_tet4(:,1),:)+VX(E_tet4(:,2),:));
%         VX_6 =0.5.*(VX(E_tet4(:,2),:)+VX(E_tet4(:,3),:));
%         VX_7 =0.5.*(VX(E_tet4(:,3),:)+VX(E_tet4(:,1),:));
%         VX_8 =0.5.*(VX(E_tet4(:,1),:)+VX(E_tet4(:,4),:));
%         VX_9 =0.5.*(VX(E_tet4(:,2),:)+VX(E_tet4(:,4),:));
%         VX_10=0.5.*(VX(E_tet4(:,3),:)+VX(E_tet4(:,4),:));
%         VX10=[VX_1_4;VX_5;VX_6;VX_7;VX_8;VX_9;VX_10];        
%         V_tet10_cell{q}=VX10(ind_uni_1,:);
%     end
% else
%     V_tet10_cell={};
% end

%% Compose output
varargout{1}=E_tet10;
varargout{2}=V_tet10;
varargout{3}=V_tet10_cell;
varargout{4}=Fb_tet10;
varargout{5}=S;
   
 
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
