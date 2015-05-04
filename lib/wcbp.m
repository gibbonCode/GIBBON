function cMap=wcbp(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cPink=[220 000 163];
cBlue=[0 0 204];
cCyan=[0 160 160];

cMap=[255*ones(1,3);... 
    120*ones(1,3);...
    cCyan;...
    cBlue;...
    cPink/1.5;...
    cPink;]./255;

[cMap]=resampleColormap(cMap,n);