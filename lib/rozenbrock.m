function [Z]=rozenbrock(X,Y)

% -----------------------------------------------------------------------
% function [Z]=rozenbrock(X,Y)
%
% "In mathematical optimization, the Rosenbrock function is a non-convex
% function used as a test problem for optimization algorithms. It is also
% known as Rosenbrock's valley or Rosenbrock's banana function. This
% function is often used to test performance of optimization algorithms.
% The global minimum is inside a long, narrow, parabolic shaped flat
% valley. To find the valley is trivial, however to converge to the global
% minimum is difficult.It is defined by Z=(1-X.^2)+100.*((Y-(X.^2)).^2)." 
% 
% From: http://en.wikipedia.org/wiki/Rosenbrock_function
%
% Kevin Moerman
% kevinmoerman@hotmail.com
% -----------------------------------------------------------------------

Z=(1-X.^2)+100.*((Y-(X.^2)).^2);
 
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
