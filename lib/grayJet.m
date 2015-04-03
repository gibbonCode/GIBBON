function cMap=grayJet(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

ns=7;
cMap=[flipud(gray(ns)); jet(ns)];

[cMap]=resampleColormap(cMap,n);