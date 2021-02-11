function M=sphereIm(varargin)

% function M=sphereIm(n,logicFlag)
%-------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------

%%
switch nargin
    case 1
        n=varargin{1};
        logicFlag=1;
    case 2
        n=varargin{1};
        logicFlag=varargin{2};
end

%%

[I,J,K]=ndgrid(-n:n,-n:n,-n:n);
M=sqrt(I.^2+J.^2+K.^2);

if logicFlag==1
    M=M<=n;
end

