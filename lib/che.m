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
cMap=[122 141 155; 0 49 79; 214 20 28; 254 229 163; ]./255; 

%% Resample colormap
[cMap]=resampleColormap(cMap,n);

end
 
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
