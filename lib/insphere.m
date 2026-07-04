function [R,Xc]=insphere(E,X)
% function [R,Xc]=insphere(E,X)
% ------------------------------------------------------------------------
% Computes the incentres Vc, and inraddi R of the inspheres for the
% tetrahedral input mesh defined by the element array E and the vertices V.
% ------------------------------------------------------------------------

TR = triangulation(E,X); %Get triangulated object
Xc = incenter(TR); %Calculate incenter coordinates
V  = tetVol(E,X); %Get element volumes

%Compute element areas
A1 = patchArea(E(:,[1 2 3]),X);
A2 = patchArea(E(:,[1 2 4]),X);
A3 = patchArea(E(:,[2 3 4]),X);
A4 = patchArea(E(:,[3 1 4]),X);
A  = A1 + A2 + A3 + A4; %Element areas

R  = 3*V./A;

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
