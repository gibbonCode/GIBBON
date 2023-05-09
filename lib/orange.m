function [c]=orange(varargin)

c=[255 123 21]/255;
switch nargin 
    %case 0
    %c=[255 123 21]/255;
    case 1
        m=varargin{1};
        c=repmat(c,m,1);
end
