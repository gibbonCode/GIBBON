function [varargout]=element2patch(varargin)

%%
% ------------------------------------------------------------------------
% This function generates FACE data for patch grapics display of elements
%
%N.B. This function assumes a node order whereby 1 2 3 4 described the top
%face and 5 6 7 8 the bottom face for hexahedral elements
% ------------------------------------------------------------------------

%% PARSE INPUT

switch nargin
    case 1
        E=varargin{1};
        C=[];
        elementType=[];
    case 2
        E=varargin{1};
        C=varargin{2};
        elementType=[];
    case 3
        E=varargin{1};
        C=varargin{2};
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
        case 4 %Linear tets
            elementType='tet4';
        case 10 %Quadratic tets
            elementType='tet10';
        case 6 %Quadratic triangles
            elementType='tri6';
        case 3
            elementType='tri3';
    end    
    disp([elementType,' elements assumed, for other elements please specify elementType']);
end

%%

switch elementType
    case 'tri3' %Linear triangles
        F=E;
        C=C;
    case 'tri6' %Quadratic triangles
        F=E(:,[1 4 2 5 3 6]);            
        C=C; 
    case 'quad4' %Linear quadrangles
        F=E;
        C=C;
    case 'quad8' %Quadratic quadrangles
        %TO DO
    case 'tet4' %Linear tets
        F=[E(:,[2 1 3]);... %face 1 2 3
            E(:,[1 2 4]);... %face 1 2 4
            E(:,[2 3 4]);... %face 2 3 4
            E(:,[3 1 4])];   %face 1 3 4
        F=fliplr(F);
        C=repmat(C,4,1); %Replicate color data
    case 'tet10' %Quadratic tets
         F=[E(:,[2 5 1 7  3 6 ]);... %face 1 2 3
            E(:,[1 5 2 9  4 8 ]);... %face 1 2 4
            E(:,[2 6 3 10 4 9 ]);... %face 2 3 4
            E(:,[3 7 1 8  4 10])];   %face 1 3 4        
        F=fliplr(F);
        C=repmat(C,4,1); %Replicate color data
    case 'hex8' %Hexahedral elements
        F =[E(:,[4 3 2 1]);... %top
            E(:,[5 6 7 8]);... %bottom
            E(:,[1 2 6 5]);... %side 1
            E(:,[3 4 8 7]);... %side 2
            E(:,[2 3 7 6]);... %front
            E(:,[5 8 4 1]);]; %back                
        C=repmat(C,6,1);
    case 'hex20' %Hexahedral elements
        %TO DO
end

%% Compose output
varargout{1}=F; 
varargout{2}=C; 

end