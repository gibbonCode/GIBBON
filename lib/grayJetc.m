function cMap=grayJetc(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[1 0 0; 1 1 0; 0 1 0; 0 1 1; ];
cMap=[flipud(gray(size(cMap,1))); flipud(cMap)];
        

[cMap]=resampleColormap(cMap,n);