%% crossProdMat
% Below is a demonstration of the features of the |crossProdMat| function

%%
clear; close all; clc;

%% Syntax
% |A=crossProdMat(a);|

%% Description 
% This function computes the so-called "cross product matrix". The output
% is a matrix A which is a skew-symmetric tensor, which allows for the
% computation of the cross product with the input vector a. In other words
% the matrix A can be used to compute c=A*b which is equivalent to
% c=cross(a,b).

%% Examples 
% 

% Example vector
a=[1 0 0]'

% Compute the cross product matrix
[A]=crossProdMat(a)

% Example vector to compute cross product with 
b=[0 1 0]'

% Use A to compute cross product axb
c=A*b

% Check "traditional" approach
c=cross(a,b)

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
