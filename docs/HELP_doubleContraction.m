%% doubleContraction
% Below is a demonstration of the features of the |doubleContraction| function

%%
clear; close all; clc;

%% Syntax
% |C=doubleContraction(A,B,subContract);|

%% Description 
% This function performs a double contraction of two tensors |A| and |B| to
% produce the output tensor |C|. The tensors are either both 3x3 second
% order tensors or |A| is a 4th order 3x3x3x3 tensor. The optional 3rd
% input |subContract| defines the indices to contract over for 4th order
% tensors e.g. |subContract=[1 2]| would contract over the first 2 indices.
% 

%% Examples 
% 

%% Double contraction of a 2nd-order and a 2nd-order tensor

%%
% Creating an example tensor. Here a tensor |A| with all zeros except for
% the 3,3 component is created. Next a tensor |B| which, through double
% contraction, can "pick out" this scalar value is created. Finally the
% double contraction is performed illustrating this process. 

%Create example tensor A
A=zeros(3,3); A(3,3)=pi; %Zeros with Azz=A33=pi
R=euler2DCM([0.25*pi 0.25*pi 0]); %Rotation matrix
A=R*A*R' %Rotated version of A

% Create 'structure' tensor example B
b=[0 0 1]*R.'; %Rotated Z vector
B=b.'*b %Rotated "structure" tensor for vector b

%%
% Double contraction using |doubleContraction|
s=doubleContraction(A,B)

%% Double contraction of a 4th-order and a 2nd-order tensor

%%
% Example for Hooke's law of linear elasticity

%Creating the 4th order linear elastic stiffness tensor
I=eye(3,3); %The 2nd order identity tensor
II1=dyadicProduct(I,I,1); %4th order base tensor 1                                                                
II3=dyadicProduct(I,I,3); %4th order base tensor 3

%Parameters for Hooke's law
mu=1; %The shear modulus
lambda=5; %The lambda lame parameter

%Compose stiffness tensor
C=lambda.*II1+2.*mu.*II3; 

%Create strain tensor
E=rand(3,3); E=E*E'

%%
% Double contraction using |doubleContraction|
subContract=[3 4];
S=doubleContraction(C,E,subContract)

%% Using symbolic variables

%%
% Example for Hooke's law of linear elasticity

%Creating the 4th order linear elastic stiffness tensor
I=eye(3,3); %The 2nd order identity tensor
II1=dyadicProduct(I,I,1); %4th order base tensor 1                                                                
II3=dyadicProduct(I,I,3); %4th order base tensor 3

%Parameters for Hooke's law
try
    syms mu lambda
catch
    mu=1; %The shear modulus
    lambda=5; %The lambda lame parameter
end
%Compose stiffness tensor
C=lambda.*II1+2.*mu.*II3;

%Display Voigt mapped version
c=voigtMap(C)

%Create strain tensor
syms E1_1 E2_2 E3_3 E1_2 E1_3 E2_3;
E=[E1_1 E1_2 E1_3;...
   E1_2 E2_2 E2_3;...
   E1_3 E2_3 E3_3]

%%
% Double contraction using |doubleContraction|

subContract=[3 4];
S=doubleContraction(C,E,subContract)

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
