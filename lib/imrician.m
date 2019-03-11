function M_rice=imrician(M,s)

% function M_rice=imrician(M,s)
% ------------------------------------------------------------------------
% IMCRICIAN Random samples from the Rice/Rician probability distribution.
%
% R ~ Rice(v, s) if R = sqrt(X^2 + Y^2), where X ~ N(v*cos(a), s^2) and
% Y ~ N(v*sin(a), s^2) are independent normal distributions (any real a).
%
% Reference: http://en.wikipedia.org/wiki/Rice_distribution
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/04/2009
% ------------------------------------------------------------------------

X_GAUSS=s.*randn(size(M))+M;
Y_GAUSS=s.*randn(size(M));
M_rice = hypot(X_GAUSS,Y_GAUSS); 
 
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
