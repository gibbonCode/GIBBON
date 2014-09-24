function [varargout]=meshBoundary(varargin)

switch nargin
    case 1
        E=varargin{1};
        [F]=element2patch(E); %Get mesh faces
    case 2
        E=varargin{1};
        elementType=varargin{2};
        [F]=element2patch(E,[],elementType); %Get mesh faces
    case 3        
        E=varargin{1};        
        elementType=varargin{2};        
        C=varargin{3};
        [F,CF]=element2patch(E,C,elementType); %Get mesh faces
    otherwise
        error('Wrong number of input arguments');
end

%Check shared faces faces
numFacesIni=size(F,1);
[F_uni,indF,IND_F_2]=uniqueIntegerRow(F);
CF=CF(indF,:);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Compose boundary set from faces that are used once
logicBoundary=F_count==1; 
Fb=F_uni(logicBoundary,:);
CFb=CF(logicBoundary,:);

varargout{1}=Fb;
varargout{2}=F_count;
varargout{3}=CFb;
