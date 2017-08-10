%% viewFourthOrderTensor
% Below is a demonstration of the features of the |viewFourthOrderTensor2D| function

%%
clear; close all; clc;

%% Syntax
% |viewFourthOrderTensor(C,numDigits,fontSizeIm,fontSize);|

%% Description 
% This function creates the 9x9 and the 6x6 (aka Voigt) array mappings for
% 3D 4th order input tensor C. The function can take up to 4 inputs: 
%%
% *  |C|, a fourth order tensor (3x3x3x3)
% *  |numDigits|, the number of decimal places to use for numerical values (default 5)
% *  |fontSizeIm|, the font size in the image (default 15)
% *  |fontSize|, the font size of the axis window (default 15)

%% Examples 

%% Viewing fourth-order stiffness tensors

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
% Visualizing the tensor using |viewFourthOrderTensor|

viewFourthOrderTensor(C); %Visualize tensor C

%% Viewing fourth-order stiffness tensors with symbolic variables

syms mu lambda; %Create symbolic Lame parameters
C=lambda.*II1+2.*mu.*II3; %Construct 4th order stiffness tensor

%%
% Visualizing the tensor using |viewFourthOrderTensor|
viewFourthOrderTensor(C); %Visualize tensor C

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
