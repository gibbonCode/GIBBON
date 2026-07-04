%% subTet
% Below is a demonstration of the features of the |subTet| function

%% Syntax
% |[Es,Vs]=subTet(E,V,splitMethod);|

%% Description 
% 
%% Examples 

clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceColor1='g';
faceColor2='r';
faceAlpha1=0.3;
faceAlpha2=1;
edgeColor=0.*ones(1,3);
edgeWidth=2;
markerSize=2;
cMap=gjet(250);

%% 
% Create a test tesselation
r=1; %Radius of tetrahedron circumsphere
[V,~]=platonic_solid(1,r);
E=[1:4];

%% Example creating sub-tetrahedrons

%Method 1
[E1,V1]=subTet(E,V,1);
 
%Method 2
[E2,V2]=subTet(E,V,2);

%% Visualization

C1=(1:size(E1,1))'; %Element colors
[F1,CF1]=element2patch(E1,C1); %Patch data for plotting

C2=(1:size(E2,1))'; %Element colors
[F2,CF2]=element2patch(E2,C2);  %Patch data for plotting

cFigure;
subplot(1,2,1); 
title('Method 1','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F1,'Vertices',V1,'FaceColor','flat','CData',CF1,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
set(gca,'FontSize',fontSize);

subplot(1,2,2); 
title('Method 2','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F2,'Vertices',V2,'FaceColor','flat','CData',CF2,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

colormap(cMap); colorbar;
view(3); grid on; axis equal; axis tight;
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
