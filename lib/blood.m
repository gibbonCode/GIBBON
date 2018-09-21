function cMap=blood(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[1 0.5 0.4; 0.9 0.3 0.27; 0.8 0.2 0.18; 0.7 0.1 0.09; 0.6 0 0; 0.5 0 0; 0.4 0 0;];

[cMap]=resampleColormap(cMap,n);