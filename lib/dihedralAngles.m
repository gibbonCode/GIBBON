function [A]=dihedralAngles(varargin)

% function [A]=dihedralAngles(E,V,elementType)
% -----------------------------------------------------------------------
% This function computes the dihedral angles A for the input elements
% defined by E, the element nodal connectivity matrix, and V, the node or
% vertex coordinate matrix. The output array A contains the same number of
% rows as E (a row for each element), and contains m columns where m is the
% number of edges for this type of element (e.g. 12 for a hexahedral element
% or 6 for a tetrahedral element). 
%
%
% 2021/03/24 Created KMM
% -----------------------------------------------------------------------

%% PARSE INPUT

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        elementType=[];
    case 3
        E=varargin{1};
        V=varargin{2};
        elementType=varargin{3};
    otherwise
        error('Wrong number of inputs');
end

numNodes=size(E,2);

if isempty(elementType) %have to assume defaults    
    switch numNodes             
        case 8 %Hexahedral elements
            elementType='hex8';            
        case 20 %Hexahedral elements
            elementType='hex20';
            numFaces=6; 
        case 4 %Linear tets
            elementType='tet4';
        case 10 %Quadratic tets
            elementType='tet10';
            numFaces=4; 
        case 6 %Quadratic triangles
            elementType='tri6';            
        case 3
            elementType='tri3';            
        case 14
            elementType='rhomdo14';            
        otherwise            
            elementType='unknown';            
    end    
end

switch elementType
    case 'tet4'
        numFaces=4;
        c  = [1 3; 3 2; 2 1;...
              1 4; 2 4; 3 4];
        
        ee = [3 2; 4 3; 1 3;...
              1 2; 4 1; 2 4];
    case 'hex8'
        numFaces=6;
        c=[1 3; 1 4; 1 5; 1 6; ...
           2 3; 2 4; 2 5; 2 6; ...
           3 5; 3 6; 4 5; 4 6; ...
           ];
        
        ee=[2 1; 4 3; 3 2; 1 4;...
            5 6; 7 8; 6 7; 8 5;...
            2 6; 5 1; 7 3; 4 8;...
            ];
end

%%

[F,~]=element2patch(E,[],elementType);
N=patchNormal(F,V);
  
FE=reshape(1:1:size(F,1),size(E,1),numFaces);
A=zeros(size(E,1),size(c,1));
for q=1:1:size(c,1)
    a=acos(dot(N(FE(:,c(q,1)),:),N(FE(:,c(q,2)),:),2));
    p1=V(E(:,ee(q,1)),:);
    p2=V(E(:,ee(q,2)),:);    
    n=vecnormalize(p2-p1); 
    x=cross(N(FE(:,c(q,1)),:),N(FE(:,c(q,2)),:),2);
    logic1=dot(n,x,2)>0;
    logicAngle=(a<pi & logic1);
    a(logicAngle)=pi-a(logicAngle);
    A(:,q)=a;
end
