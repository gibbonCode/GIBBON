function [y]=smoothSquareWave(varargin)

% function [y]=smoothSquareWave(t,d)
% ------------------------------------------------------------------------
% 
% ------------------------------------------------------------------------

%%

switch nargin
    case 1
        t=varargin{1};
        d=0.1;
    case 2
        t=varargin{1};
        d=varargin{2};
end

%%
y=atan(sin(t)/d)./atan(1/d);

