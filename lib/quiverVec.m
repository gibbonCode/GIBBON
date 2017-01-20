function [varargout]=quiverVec(varargin)


%% Parse input

switch nargin
    case 2
        P=varargin{1};
        V=varargin{2};
        vecSize=[];
        colorOpt=[];
        edgeColorOpt='none';
    case 3
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorOpt=[];
        edgeColorOpt='none';
    case 4
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorOpt=varargin{4};        
        edgeColorOpt='none';
    case 5
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorOpt=varargin{4};
        edgeColorOpt=varargin{5};
end


if isempty(vecSize)
    [F,P,C]=quiver3Dpatch(P(:,1),P(:,2),P(:,3),V(:,1),V(:,2),V(:,3),[],[]);
else
    [F,P,C]=quiver3Dpatch(P(:,1),P(:,2),P(:,3),V(:,1),V(:,2),V(:,3),[],[vecSize(1) vecSize(2)]);
end

if isempty(colorOpt) %If empty use colormapping
    h=gpatch(F,P,C,edgeColorOpt,1);
else %else use specified which could be 'k'
    h=gpatch(F,P,colorOpt,edgeColorOpt,1);
end

varargout{1}=h;