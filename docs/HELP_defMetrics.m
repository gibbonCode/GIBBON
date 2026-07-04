%% defMetrics
% Below is a demonstration of the features of the |defMetrics| function

%%
clear; close all; clc;

%% Syntax
% |[D_out]=defMetrics(F_cell,strain_type);|

%% Description 
% Computes various deformation metrics using the input cell array F_cell
% which contains deformation gradient tensors.

%% Examples 
% 

F_uni = [1.5 0.0 0.0;...
         0.0 1.0 0.0;...
         0.0 0.0 1.0];

F_shear=[1.0 0.5 0.0;...
         0.0 1.0 0.0;...
         0.0 0.0 1.0];

F_hydro=[1.5 0.0 0.0;...
         0.0 1.0 0.0;...
         0.0 0.0 1.0];

F_rand=1+0.25*(rand(3,3)-0.5);

F_cell={F_uni,F_shear,F_hydro,F_rand};

% 1 = Biot (linear) strain tensor
% 2 = Hencky (logarithmic/natural) strain tensor
% 3 = Green-Lagrange strain tensor
strain_type=2; 
[D_out]=defMetrics(F_cell,strain_type);

D_out

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
