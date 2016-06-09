function cMap=kwJet(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[0 0 0; gjet(4); 1 1 1];

[cMap]=resampleColormap(cMap,n);