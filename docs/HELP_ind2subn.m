%% ind2subn
% Below is a demonstration of the features of the |ind2subn| function

%%
clear; close all; clc;

%% Syntax
% |[A] = ind2subn(siz,ind);|

%% Description 
% This function is similar to MATLAB's ind2sub function. However the output
% is collected in an array that is numel(ind) x numel(siz). This function
% constructs an array where the columns represent the subscripts indices
% for the input linear indices (ind). 
%
%   See also: |ind2sub|, |sub2ind|

%% Examples 
% 

%%
% 2D indices
ind=[1 2 6 12]
siz=[6 6]
[IJ] = ind2subn(siz,ind)

%%
% 3D indices
ind=[1 2 6 12 666]
siz=[11 11 11]
[IJK] = ind2subn(siz,ind)

%%
% 4D indices
ind=[1 2 6 12 666 6666]
siz=[11 11 11 11]
[IJKL] = ind2subn(siz,ind)

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
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
