%% ellipseFit
% Below is a demonstration of the features of the |ellipseFit| function

%%
clear; close all; clc;

%% Syntax
% |[A] = ellipseFit(V,optMethod,numSample);|

%% Description 
% 

%% Examples 
% 

n=25;
t=linspace(0,2*pi,n+1); 
t=t(1:end-1);

x=2*cos(t);
y=3*sin(t);
V=[x(:) y(:) zeros(size(x(:)))]+5*randn(1,3);
V=V+0.1*randn(size(V));

%% Fit ellipse

% Fit ellipse
[A] = ellipseFit(V);

% Get coordinates for plotting
[VF]=ellipseCoord(A,linspace(0,2*pi,500));

%%
cFigure;hold on;
plotV(VF,'b-','LineWidth',2);
plotV(V,'r.','MarkerSize',35);
axis equal; box on;
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
