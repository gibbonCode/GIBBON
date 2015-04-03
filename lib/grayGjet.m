function cMap=grayGjet(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=gjet(4);
cMap=[flipud(gray(size(cMap,1))); cMap];
        

[cMap]=resampleColormap(cMap,n);