%% cellEig
% Below is a demonstration of the features of the |cellEig| function

%% Syntax
% |[V,D]=cellEig(C);|

%% Description 
% Computes eigenvalues and eigenvectors for each matrix contained in the
% cell array C, i.e. [v,d]=eig(c) is executed for each cell entry. The
% output is two cell arrays, i.e. the cell V containing the eigenvectors
% and the cell D containing the eigenvalues. 

%% Examples

%%
clear; close all; clc;

%% Example: Calculating eigenvalues for matrices contained in cells
% Creating example cell containing two matrices

M1=rand(3,3);
M1=M1*M1';
M2=rand(5,5);
M2=M2*M2';

C={M1,M2};
[V,D]=cellEig(C);

%%
% Contained in the output cells are the eigenvectors and eigenvalues of
% each of the matrices e.g. for the first
v1=V{1}
d1=D{1}

%%
% and the second entry
v2=V{2}
d2=D{2}

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
