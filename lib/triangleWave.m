function s=triangleWave(varargin)

switch nargin
    case 1
        t=varargin{1};
        p=2*pi;
    case 2
        t=varargin{1};
        p=varargin{2};
end
s=sawtoothWave(t+p/4,p);
s=abs(2*s)-1;


