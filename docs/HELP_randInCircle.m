%% randInCircle
% Below is a demonstration of the features of the |randInCircle| function

%%
clear; close all; clc;

%% Syntax
% |a=randInCircle(siz);|

%% Description 
% This function generates a matrix or array of uniformly distributed random
% points inside a circle with radius R. 

%%
% Plot settings for examples
markerSize=35; 
lineWidth=3; 
fontSize=35;

%% Example
% 

n=1000; %Number of points
R=2; %Radius of the circle
V=randInCircle(n,R); %Uniformly sampled points inside the circle

%% 
% Visualization

%Create circle coordinates to visualize circle curve
t=linspace(0,2*pi,100)';
vc=R.*[cos(t) sin(t)];

cFigure; hold on; 
hp1=plotV(V,'k.','MarkerSize',markerSize);
hp2=plotV(vc,'b-','LineWidth',lineWidth);
axis tight; axis equal; box on; grid on; 
view(2);
set(gca,'FontSize',fontSize);
legend([hp1 hp2],{'Uniformly distributed circle interior points','Circle boundary'},'Location','NorthOutside');
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
