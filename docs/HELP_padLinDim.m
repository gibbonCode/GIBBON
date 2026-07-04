%% padLinDim
% Below is a demonstration of the features of the |padLinDim| function

%%
clear; close all; clc;

%% Syntax
% |[XP,indOriginal]=padLinDim(X,numPad,padDim);|

%% Description 
% Pads an array |X| allong dimension |padDim| with |numPad| entries which
% are linearly extrapolated allong that direction.

%% Examples 
% 

%%
% Plot settings
fontSize=35; 

%% Linearly padding a vector

%%
% Create example vector
t=linspace(0,2*pi,50);
y=10*sin(t);

%%
% Linearly pad vector

numPad=5; %Number of elements to pad
padDim=2; %Dimension/direction to pad in 
[yp,indOriginal]=padLinDim(y,numPad,padDim);
[tp]=padLinDim(t,5,2);

%%
cFigure; hold on;
hp1=plot(tp,yp,'k-','LineWidth',3);
hp2=plot(t,y,'g-','LineWidth',5);
legend([hp2 hp1],{'Original','Padded'});
set(gca,'FontSize',fontSize);
axis tight; axis square; box on; grid on; 
drawnow; 

%% Linearly padding a matrix

%%
% Create example matrix
Y=y(ones(1,25),:);

%%
% Linearly pad matrix

numPad=15; %Number of elements to pad
padDim=2; %Dimension/direction to pad in 
[Yp,indOriginal]=padLinDim(Y,numPad,padDim);

%%

cFigure; hold on;
hp1=surf(Yp);
set(hp1,'FaceColor','w');
set(gca,'FontSize',fontSize);
axisGeom; camlight headlight;
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
