function S=smoothHeaviside(varargin)

switch nargin
    case 1
        x=varargin{1};
        k=6;
        r=0;
    case 2
        x=varargin{1};
        k=varargin{2};
        r=0;
    case 3
        x=varargin{1};
        k=varargin{2};
        r=varargin{3};
end
k=k*2;
S=exp(2*k*(x-r))./(1+exp(2*k*(x-r)));

