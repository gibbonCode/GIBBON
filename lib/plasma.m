function [cMap]=plasma(varargin)

% function [cMap]=plasma(n)
% ------------------------------------------------------------------------
% Creates the colormap data for the plasma colormap. 
%
% This colormap is a MATLAB implementation of the matplotlib colormap by
% Nathaniel J. Smith, Stefan van der Walt, and Eric Firing. 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/04/10
%------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%

[cMap]=matplotlibColormap('plasma',n);

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
