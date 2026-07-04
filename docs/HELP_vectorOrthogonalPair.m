%% vectorOrthogonalPair
% Below is a demonstration of the features of the |vectorOrthogonalPair| function

%%
clear; close all; clc;

%% Syntax
% |[a,d]=vectorOrthogonalPair(f);|

%% Description 
% Based on the input vector f this function generates the normalized output
% vectors a and d which are orthogonal to eachother and to f. 
% This function can be useful for describing local element axis systems
% based e.g. in FEBio (whereby the input vector defaults to e3). 

%% Examples

%% Example 1: Creating a triplet of mutually orthogonal vectors 

%%
% Creating example vectors

P=eye(3,3); %Vector origins
V=euler2DCM([0.25*pi 0.25*pi 0.25*pi]); %Input Vectors, rotated directions
% V=eye(3,3); %Input Vectors, x,y,z axes

%%
% Compute mutually orthogonal sets using |vectorOrthogonalPair|
[a,d]=vectorOrthogonalPair(V);

%%
% Visualize the sets

cFigure; 
title('Visualized vectors sets, input vectors=transparent, orthogonal=non-transparent');
quiverVec(P,V,2,'b','none',1,0.2);
quiverVec(P,a,1,'r','none',1,1);
quiverVec(P,d,1,'g','none',1,1);
% quiverVec(P,d,1,'k','none',1,1);
axisGeom; 
camlight headlight; 
drawnow; 

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
