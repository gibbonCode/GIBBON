function [cMap]=che(varargin)

% function [cMap]=che(n)
% ------------------------------------------------------------------------
% Creates the colormap data for n levels for the che colormap. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/04/14
%------------------------------------------------------------------------

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%% Define colors
cMap=[0 49 79; 214 20 28; 122 141 155; 254 229 163; ]./255; 

%% Resample colormap
[cMap]=resampleColormap(cMap,n);

end