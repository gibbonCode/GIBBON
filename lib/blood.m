function cMap=blood(varargin)

% function cMap=blood(n)
% -----------------------------------------------------------------------
% This function creates the colormap data |cMap| for the blood colormap
% using |n| levels. If |n| is not provided the default used is 250. 
%
% 2019/06/27 Updated description
% -----------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%
cMap=[1   0.6 0.48;...
      1   0.5 0.4 ;...
      0.9 0.3 0.27;...      
      0.7 0.1 0.09;...
      0.6 0   0;...
      0.3 0   0;];
  
[cMap]=resampleColormap(cMap,n);

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
