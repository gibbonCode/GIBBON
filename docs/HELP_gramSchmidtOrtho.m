%% gramSchmidtOrtho
% Below is a demonstration of the features of the |gramSchmidtOrtho| function

%%
clear; close all; clc;

%% Syntax
% |[E]=gramSchmidtOrtho(Q);|

%% Description 
% This function uses the Gram-Schmidt method to convert the input basis Q
% to an orthonormal basis E. 

%% Examples 
% 

%%
% Create example non-orthonormal base vector set

rz=euler2DCM([0 0 0.25*pi]); %Rotation around z
rx=euler2DCM([-0.25*pi 0 0]); %Rotation around x
Q=eye(3,3); %Identify matrix (orthonormal)
Q(:,1)=rz*Q(:,1); % Rotate the 1st axis towards 2nd
Q(:,3)=rx*Q(:,3); % Rotate the 3rd towards 1-2-plane

%Rotate all
R=euler2DCM([0.25*pi 0.25*pi 0.25*pi]); 
Q=R*Q;

%Scale all (so not normal)
Q=Q*2;

%%
% Use Gram-Schmidt method to obtain an orthonormal basis

[E]=gramSchmidtOrtho(Q)

%%
% Visualization

cFigure; 

subplot(1,2,1); hold on; 
title('A non-orthonormal set');
hp1=quiverTriad(zeros(1,3),Q,2);
axisGeom; camlight headlight; 

subplot(1,2,2); hold on; 
title('An orthonormal set');
hp2=quiverTriad(zeros(1,3),E,1);
axisGeom; camlight headlight; 

gdrawnow;

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
