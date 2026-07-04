function [IND_F,IND_V]=patchIND(F,V,formOpt)

% function [IND_F,IND_V]=patchIND(F,V,formOpt)
% ------------------------------------------------------------------------
%
% This function computes the neighbouring faces (read elements in the case
% of >3D tesselations) and vertices for the input tesselation specified by
% F and V. The function generates matrices whereby row entries correspond
% to vertices and column entries correspond to neighbouring (according to
% the connectivity in the tesselation specified) vertices in the case IND_V
% and faces in the case of IND_F. 
% If formOpt==2 the output arrays IND_F and IND_V are sparse arrays of size
% [size(V,1),size(V,1)] and [size(V,1),size(F,1)] respectively. However if
% formOpt==1 (DEFAULT) or not specified then the output array is instead a
% full array whose size in the column direction is smaller and depends on
% the maximum number of vertex and face neighbours encountered. E.g. in a
% triangulated mesh where each vertex is connected to 6 vertices the array
% IND_V is thus [size(V,1),6] in size. If some locations have less than
% this maximum number of vertices zeros are used to fill up the array. The
% same line of arguments holds for the IND_F array where instead entries
% reflect neighbouring faces. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/01/26 %added option of outputting as sparse arrays
%------------------------------------------------------------------------

%%

%Check if formOpt is provided
if nargin==2 %DEFAULT
    formOpt=2;  %Output is cropped and full  
end

switch formOpt
    case 1
        sparseOpt=1;
    case 2
        sparseOpt=0;
end

[IND_F,IND_V]=tesIND(F,V,sparseOpt);
 
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
