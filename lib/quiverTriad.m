function [varargout]=quiverTriad(varargin)


%% Parse input

switch nargin
    case 2
        V=varargin{1};
        R=varargin{2};
        vecSize=1;
        colorOpt=[];
    case 3
        V=varargin{1};
        R=varargin{2};
        vecSize=varargin{3};
        colorOpt=[];
    case 4
        V=varargin{1};
        R=varargin{2};
        vecSize=varargin{3};
        colorOpt=varargin{4};        
end

V=V(ones(1,3),:);

[F,V,C]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),R(:,1),R(:,2),R(:,3),eye(3,3),[vecSize vecSize]);

if isempty(colorOpt) %If empty use colormapping
    h=gpatch(F,V,C,'none',1);
else %else use specified which could be 'k'
    h=gpatch(F,V,colorOpt,'none',1);
end

varargout{1}=h;


