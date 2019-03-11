function [V,D]=cellEig(C)

% function [V,D]=cellEig(C)
% ------------------------------------------------------------------------
% Computes eigenvalues and eigenvectors for each matrix contained in the
% cell array C, i.e. [v,d]=eig(c) is executed for each cell entry. The
% output is two cell arrays, i.e. the cell V containing the eigenvectors
% and the cell D containing the eigenvalues.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 28/11/2013
% 2016/09/09 Updated documentation
%------------------------------------------------------------------------

[V,D]=cellfun(@eig,C,'UniformOutput',0);
 
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
