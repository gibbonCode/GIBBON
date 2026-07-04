%% fourthOrderCell
% Below is a demonstration of the features of the |fourthOrderCell| function

%%
clear; close all; clc;

%% Syntax
% |[CM]=fourthOrderCell(C);|

%% Description 
% Converts a 3x3x3x3 4th order tensor into a cell 3x3 cell array featuring
% 3x3 matrix entries. 

%% Examples 
% 

%%
% Creating the stiffness tensor for Hooke's law of linear elasticity

%Constructing 4th order base tensor set
I=eye(3,3); %The 2nd order identity tensor
II1=dyadicProduct(I,I,1); %4th order base tensor 1                                                                
II3=dyadicProduct(I,I,3); %4th order base tensor 3

%Parameters for Hooke's law
mu=1; %The shear modulus
lambda=5; %The lambda lame parameter
C=lambda.*II1+2.*mu.*II3; %Construct 4th order stiffness tensor

%%

[CM]=fourthOrderCell(C)

CM{1,1}

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
