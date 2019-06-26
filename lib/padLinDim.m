function [varargout]=padLinDim(varargin)

% function [XP,indOriginal]=padLinDim(X,numPad,padDim)
% -----------------------------------------------------------------------
% Pad arrays (including vectors) allong dimension padDim with numPad
% entries which are linearly extrapolated allong that direction.
%
% -----------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        X=varargin{1};
        numPad=1;
        padDim=1;
        dirOpt='both';
    case 2
        X=varargin{1};
        numPad=varargin{2};
        padDim=1;
        dirOpt='both';
    case 3
        X=varargin{1};
        numPad=varargin{2};
        padDim=varargin{3};
        dirOpt='both';
    case 4
        X=varargin{1};
        numPad=varargin{2};
        padDim=varargin{3};
        dirOpt=varargin{4};
end

switch dirOpt
    case 'both'
        [XP,indOriginal]=padLinDimStart(X,numPad,padDim);
        XP=flip(XP,padDim);
        XP=padLinDimStart(XP,numPad,padDim);
        XP=flip(XP,padDim);
    case 'start'
        [XP,indOriginal]=padLinDimStart(X,numPad,padDim);
    case 'end'
        X=flip(X,padDim);
        [XP,indOriginal]=padLinDimStart(X,numPad,padDim);
        XP=flip(XP,padDim);
end

varargout{1}=XP;
varargout{2}=indOriginal;

end

%%

function [varargout]=padLinDimStart(X,numPad,padDim)
numDims=ndims(X);
siz=size(X);

Q=1:1:numDims;
Q=Q(Q~=padDim);
n=prod(siz(Q));
A=zeros(n,numDims);
A(:,padDim)=1;
for q=Q(:)'
    s=siz(q);
    A(:,q)=repmat(1:1:s,1,n/s);
end
B=A;
B(:,padDim)=2;

indDiff_1=sub2indn(siz,A);
indDiff_2=sub2indn(siz,B);

sizD=siz;
sizD(padDim)=1;

X1=reshape(X(indDiff_1),sizD);
X2=reshape(X(indDiff_2),sizD);
dX=X2-X1;

sizn=siz;
sizn(padDim)=sizn(padDim)+numPad;

IND=1:1:prod(sizn);
A=ind2subn(sizn,IND);
logicOld=(A(:,padDim)-numPad)>0;
indOriginal=IND(logicOld);

XP=zeros(sizn);
XP(indOriginal)=X;

A=zeros(n,numDims);
for qq=Q(:)'
    s=sizn(qq);
    A(:,qq)=repmat(1:1:s,1,n/s);
end

for q=1:1:numPad
    A(:,padDim)=q;
    indSet=sub2indn(sizn,A);
    XP(indSet)=((numPad-q)+1);
    XP(indSet)=X1-(((numPad-q)+1).*dX);
end

varargout{1}=XP;
varargout{2}=indOriginal;
end
