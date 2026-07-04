%% ind2logic
% Below is a demonstration of the features of the |ind2logic| function

%%
clear; close all; clc;

%% Syntax
% |[L]=ind2logic(ind,siz);|

%% Description 
% This function converts the linear indicies ind to the logic array L. 

%% Examples 
% 

%% Converting linear indices without specifying size
ind=[1 5 6]
[L]=ind2logic(ind)

%% Converting linear indices with size specification
ind=[1 5 6]
siz=[3,3]
[L]=ind2logic(ind,siz)

%% Converting linear indices for higher order arrays

ind=[1 5 6 26]
siz=[3,3,3]
[L]=ind2logic(ind,siz)

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
