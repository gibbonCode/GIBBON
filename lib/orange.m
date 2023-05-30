function [c]=orange(varargin)

% function [c]=orange(m)
% ------------------------------------------------------------------------
% This function outputs the rgb color values (between 0 and 1) for an
% orange color, i.e. [255 123 21]/255
% If the optional output m is provided it is used to copy the color m
% times, such than an mX3 array is returned. 
% ------------------------------------------------------------------------

%%

%Orange color
c=[255 123 21]/255;

%Replicate if needed
switch nargin 
    %case 0
    %c=[255 123 21]/255;
    case 1
        m=varargin{1};
        c=repmat(c,m,1);
end

