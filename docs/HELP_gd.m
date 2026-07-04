%% gd
% Below is a demonstration of the features of the |gd| function

%%
clear; close all; clc;

%% Syntax
% |y=gd(x);|

%% Description 
% The output of is function |y| is the Gudermannian of the input |x|. The
% Gudermannian is defined as |y=asin(tanh(x))|.  
% See also |gd|.

%% Examples 
% 

%%
% Plot setttings
lineWidth=3; 
fontSize=25;

%% Example 1: Computing the inverse Gudermannian
x=linspace(-3*pi,3*pi,500); %A range of data
y=gd(x); %inverse Gudermannian

%%
% Visualize

cFigure; 
plot(x,y,'b-','LineWidth',lineWidth);
axis square; axis tight; grid on; box on;
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
