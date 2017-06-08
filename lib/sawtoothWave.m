function s=sawtoothWave(varargin)

switch nargin
    case 1
        t=varargin{1};
        p=2*pi;
    case 2
        t=varargin{1};
        p=varargin{2};
end

s=2*(t/p-floor(1/2+t/p));
