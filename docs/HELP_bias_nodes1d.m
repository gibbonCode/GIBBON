%% bias_nodes1d
% Below is a demonstration of the features of the |bias_nodes1d| function

%% Syntax
% |[xb]=bias_nodes1d(x,f_bias);|

%% Description 
% The |bias_nodes1d| function is able to adjust the point (or node) spacing
% of a curve based on a bias factor and biasing scheme.

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
markerSize1=25;

%% Example: Biasing node spacing allong a curve

f_bias=1.8; %Bias factor
n=15; %Number of steps
x=linspace(0,10,n);
[xb]=bias_nodes1d(x,f_bias);

%%
% Plotting results
hf1=cFigure;
title('Biased nodes','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plot(x,ones(size(x)),'k+','MarkerSize',markerSize1);
plot(xb,ones(size(x)),'r.','MarkerSize',markerSize1);

axis equal; axis tight;  grid on;  
set(gca,'FontSize',fontSize);
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
