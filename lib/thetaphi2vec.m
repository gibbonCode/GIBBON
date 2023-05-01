function [r,c]=thetaphi2vec(THETA,PHI)

%--------------------------------------------------------------------------
% function [vx,vy]=thetaphi2vec(THETA,PHI)
% 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/02/2010
%-------------------------------------------------------------------------- 


%%

r = [cos(THETA).*cos(PHI)  -sin(THETA) cos(THETA).*sin(PHI)];
c = [sin(THETA).*cos(PHI)  cos(THETA)  sin(THETA).*sin(PHI)];

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
