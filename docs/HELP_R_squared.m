%% R_squared
% Below is a demonstration of the features of the |R_squared| function

%%
clear; close all; clc;

%% Syntax
% |[R_sq]=R_squared(y,yf);|

%% Description 
% This function calculates the coefficient of determination (R-Squared) for
% the (e.g. measurement) data |y| and the data |yf| (e.g. a model fit). 
% NaN entries in either inputs are ignored. 

%% Examples 
% 

%%
% Plot settings
lineWidth=3;
markerSize=50;
fontSize=25; 

%%
% Example data
t=linspace(0,2*pi,25);
y=sin(t)+0.1*randn(size(t));

%%
% Example fit
yf=sin(t); 

%% 
% Compute R-squared
[R_sq]=R_squared(y,yf)

%%
% Visualize 

cFigure; hold on;
title(['R^2=',sprintf('%.4f',R_sq)],'Interpreter','Tex');
plot(t,y,'k.','MarkerSize',markerSize); 
plot(t,yf,'g-','LineWidth',lineWidth); 
set(gca,'FontSize',fontSize);
axis tight; axis square; box on; grid on;
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
