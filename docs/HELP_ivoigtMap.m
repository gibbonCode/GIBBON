%% ivoigtMap
% Below is a demonstration of the features of the |ivoigtMap| function

%%
clear; close all; clc;

%% Syntax
% |c=ivoigtMap(cVoigt);|

%% Description 
% This function reverses Voigt mapped tensors from their Voigt array for
% back to their tensorial form. 

%% Examples 
% 

%% Example 1: The Voigt mapping of the elasticity tensor for Hooke's law with Lame parameters

%%

%Constructing 4th order base tensor set
I=eye(3,3); %The 2nd order identity tensor
II1=dyadicProduct(I,I,1); %4th order base tensor 1                                                                
II3=dyadicProduct(I,I,3); %4th order base tensor 3

%Lame parameters for Hooke's law
mu=2; %The shear modulus
lambda=3; %The lambda lame parameter

%Construct 4th order stiffness tensor
C=lambda.*II1+2.*mu.*II3; 

Cv=voigtMap(C) % Voigt map version

Cr=ivoigtMap(Cv) 

%%

s=rand(3,3);
S=s*s'

Sv=voigtMap(S)

Sr=ivoigtMap(Sv)
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
