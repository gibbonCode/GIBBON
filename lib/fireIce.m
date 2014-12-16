function cmap=fireIce(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cmap=[1 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 1];

[cmap]=resampleColormap(cmap,n);