function [V]=ellipseCoord3(e,t)

% function [V]=ellipseCoord3(e,t)
% ------------------------------------------------------------------------
% Calculates ellipse coordinates for the angles in t based on the input
% structure e which contains the following fields: 
% radii, a 2x1 array
% axes, a 3x3 rotation matrix
% centre the ellipse centre coordinates% 
% 
% ------------------------------------------------------------------------

%%

V=[e.radii(1).*cos(t(:)) e.radii(2).*sin(t(:)) zeros(numel(t),1)];
V=(e.axes*V')';
V=V+e.centre(ones(numel(t),1),:);
 
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
