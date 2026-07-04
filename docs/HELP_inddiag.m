%% inddiag
% Below is a demonstration of the features of the |inddiag| function

%%
clear; close all; clc;

%% Syntax
% |[ind]=inddiag(A);|

%% Description 
% Returns the indices of the diagonal elements of A. If A is not a 2D array
% an error is return. If size(A,1)~=size(A,2) the function still returns
% the valid indices for when i=j in terms of A_ij. 

%% Examples 
%

%%
% 

[ind]=inddiag(rand(1,1))

%%
%
[ind]=inddiag(rand(2,2))

%%
%
[ind]=inddiag(rand(5,5))

%%
%
[ind]=inddiag(rand(5,3))

%%
%
[ind]=inddiag(rand(3,5))

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
