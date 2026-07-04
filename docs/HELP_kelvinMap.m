%% kelvinMap
% Below is a demonstration of the features of the |kelvinMap| function

%%
clear; close all; clc;

%% Syntax
% |cKelvin=kelvinMap(c);|

%% Description 
% This function creates the Kelvin mapped (6x6) tensor form of the 4-th
% order (3x3x3x3) input tensor. 

%% Examples 
% 

%% Example 1: The Kelvin mapping of the elasticity tensor for Hooke's law with Lame parameters

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

%%
% Derive Kelvin mapped tensor
Ck=kelvinMap(C) 

%%
% The Kelvin mapping for eigenstiffness determination
[V,D]=eig(Ck)

%%
% Using symbolic variables

try
    %Lame parameters for Hooke's law
    syms mu lambda; %Create symbolic parameters
    
    %Construct 4th order stiffness tensor
    C=lambda.*II1+2.*mu.*II3; 
    
    %%
    % Derive Kelvin mapped tensor
    
    Ck=kelvinMap(C)
    
    %%
    % The Kelvin mapping for eigenstiffness determination
    
    [V,D]=eig(Ck)

end

%% Example 2: The Kelvin mapping of the elasticity tensor for Hooke's law with bulk/shear modulus

% Using mu and bulk modulus parameters
mu=2; %The shear modulus
lambda=3; %The lambda lame parameter
k=lambda+2/3*mu; %Bulk modulus

%Construct 4th order stiffness tensor
C=(k-2/3*mu).*II1+2.*mu.*II3; %Construct 4th order stiffness tensor

%%
% Derive Kelvin mapped tensor

Ck=kelvinMap(C) 

%%
% The Kelvin mapping for eigenstiffness determination
[V,D]=eig(Ck)

%% 
% Using symbolic variables

try
    syms mu k; %Create symbolic parameters
    
    %Construct 4th order stiffness tensor
    C=(k-2/3*mu).*II1+2.*mu.*II3; %Construct 4th order stiffness tensor
    
    %%
    % Derive Kelvin mapped tensor
    
    Ck=kelvinMap(C)
    
    %%
    % The Kelvin mapping for eigenstiffness determination
    
    [V,D]=eig(Ck)

end

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
