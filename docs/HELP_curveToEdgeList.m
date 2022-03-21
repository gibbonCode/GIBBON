%% curveToEdgeList
% Below is a demonstration of the features of the |curveToEdgeList| function

%%
clear; close all; clc;

%% Syntax
% |[E]=curveToEdgeList(N);|

%% Description 
% UNDOCUMENTED 
%% Examples 
% 

%%
% Create example curve
t=linspace(0,2*pi,10)';
t=t(1:end-1);

V=[cos(t) sin(t) zeros(size(t))];

%%
% Get curve edges

[E]=curveToEdgeList(V)

%%
% Visualization

cFigure; 
subplot(1,2,1); hold on;
title('Curved plotted using plot command');
plotV(V,'r-','LineWidth',2,'MarkerSize',25);
axis tight; axis equal; grid on; box on; view(2); 

subplot(1,2,2); hold on;
title('Curved plotted as edges using patch command');
gpatch(E,V,'none','g',1,2);
axis tight; axis equal; grid on; box on; view(2); 
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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
