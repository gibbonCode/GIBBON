function [cMap]=ukr(varargin)

% function [cMap]=ukr(n)
% ------------------------------------------------------------------------
% Creates the colormap data for the ukr colormap. 
% Based on: https://en.m.wikipedia.org/wiki/Flag_of_Ukraine
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log: 
% 2022/02/25 Created
%------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%

cMap=resampleColormap([255 215 0; 0 87 183; ]./255,n);

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
