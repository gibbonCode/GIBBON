function cMap=bloodBone(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[1 0.95 0.8; 0.6 0 0];
[cMap]=resampleColormap(cMap,250);