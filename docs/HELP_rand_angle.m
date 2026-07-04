%% rand_angle
% Below is a demonstration of the features of the |rand_angle| function

%%
clear; close all; clc;

%% Syntax
% |a=rand_angle(siz);|

%% Description 
% This function generates a matrix or array of uniformly distributed random
% angles (of the size siz) in the range [0 2*pi]. The function operates by
% first creating angles in the range 0-pi/2. The angles are next converted
% to unit vectors in the unit circle. The X and Y components of these
% vectors are next independantly and randomly negated. This negating or
% flipping operation causes the vectors to uniformly but randomly span the
% entire cirle. Nexts the vectors are converted to an angle in the range [0
% 2*pi].  

%%
% Plot settings for examples
markerSize=25; 
lineWidth=3; 
fontSize=35;

%% Examples 
% 

siz=[500,1]; %Size of desired output
a=rand_angle(siz); % The random angles

%%

cFigure; hold on; 
xlabel('Angles'); ylabel('Count');
histogram(a,linspace(0,2*pi,6));
set(gca,'FontSize',fontSize);
drawnow; 

%% 
% Visualization

%Create circle coordinates to visualize circle curve
t=linspace(0,2*pi,100)';
vc=[cos(t) sin(t)];
N=[cos(a) sin(a)];

cFigure; hold on; 
hp1=plotV(vc,'b-','LineWidth',lineWidth);
hp2=plotV(N,'k.','MarkerSize',markerSize);
axis tight; axis equal; box on; grid on; 
view(2);
set(gca,'FontSize',fontSize);
legend([hp1 hp2],{'Circle boundary','Uniformly distributed points on circle'},'Location','NorthOutside');
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
