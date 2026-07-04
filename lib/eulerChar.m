function [eulerVal]=eulerChar(F,V)
% function [eulerVal]=eulerChar(F,V)
% ------------------------------------------------------------------------
%
% 2023/06/13 KMM: Added the use of patchCleanUnused so that unused points
% are not counted 
% ------------------------------------------------------------------------

%%

[F,V]=patchCleanUnused(F,V); %Remove unused vertices from the list
E=patchEdges(F,1); %Get the edges
eulerVal=size(V,1)-size(E,1)+size(F,1); %Compute Euler characteristic
 
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
