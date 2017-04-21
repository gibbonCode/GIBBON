function cMap=gjet(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[51  105 232;... %Blue
      0   153 37;... %Green
      238 178 17;... %Yelow
      213 15  37;... %red
    ];
cMap=cMap./255;

[cMap]=resampleColormap(cMap,n);