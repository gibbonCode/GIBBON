%% rhombicDodecahedron
% Below is a demonstration of the features of the |rhombicDodecahedron| function

%%
clear; close all; clc;

%% 
% Plot settings
figColor='w'; figColorDef='white';
fontSize=20;
faceAlpha1=0.8;
edgeColor=0.6*ones(1,3);
lineWidth1=2;
markerSize=25;

%% Creating a patch model of a rhombic dodecahedron

r=sqrt(2)/2; %Radii, the chosen level results in X,Y width of 1

[F,V]=rhombicDodecahedron(r);

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);  
title('A rhombic dodecahedron','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

hp=patch('Faces',F,'Vertices',V,'FaceColor','g','FaceAlpha',faceAlpha1,'EdgeColor',edgeColor,'lineWidth',lineWidth1,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','k','MarkerEdgeColor','none');
colormap(hsv);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; grid on; view(-10,25);
camlight('headlight'); lighting flat; 

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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
