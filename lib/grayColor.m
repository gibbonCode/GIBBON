function [C]=grayColor(varargin)

% function [C]=grayColor(colorLevels)

%%
switch nargin
    case 0
        colorLevels=0.5;
    case 1
        colorLevels=varargin{1};
end

%%
colorLevels=colorLevels(:); %Force as column
C=colorLevels(:,ones(1,3));