function [V]=ellipseCoord(A,t)

% function [V]=ellipseCoord(A,t)
% ------------------------------------------------------------------------
% Calculates ellipse coordiantes for the angles in t based on the vector A
% which defines the centre coordinates, the radii and the angle
% respectively. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2013/24/09
%------------------------------------------------------------------------

%%
x0=A(1);
y0=A(2);
x=A(3).*cos(t);
y=A(4).*sin(t);
V=[x(:) y(:) zeros(size(x(:)))];
[R,~]=euler2DCM([0 0 -A(5)]);
V=(R*V')';
V(:,1)=V(:,1)+x0;
V(:,2)=V(:,2)+y0;
 
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
