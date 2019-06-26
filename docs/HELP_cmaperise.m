%% cmaperise
% Below is a demonstration of the features of the |cmaperise| function

%%
clear; close all; clc;

%% Syntax
% |Cmapped=cmaperise(C,cmap,clim);|

%% Description 
% This function creates RGB colors for the data C using the colormap cmap
% and the limits clim. 

%% Examples 
% 

[X,Y,Z]=peaks(50);
[F,V,C]=surf2patch(X,Y,Z,Z);

clim=[min(C(:)) max(C(:))];

%RGB color set with map 1
cmap=viridis(250);
C_rgb_1=cmaperise(C,cmap,clim);

%RGB color set with map 2
cmap=gray(250);
C_rgb_2=cmaperise(C,cmap,clim);

%RGB color set with map 3
cmap=gjet(250);
C_rgb_3=cmaperise(C,cmap,clim);


%% Visualize

% Offset coordinates so plots are side by side
V2=V;
V2(:,1)=V2(:,1)-min(V2(:,1))+max(V(:,1));
V3=V2;
V3(:,1)=V3(:,1)-min(V3(:,1))+max(V2(:,1));

cFigure;

subplot(1,2,1);
title('Colormapped')
hp=gpatch(F,V,C);
legend(hp,'Colormapped data');
axisGeom;
camlight headlight; colormap viridis; colorbar;

subplot(1,2,2);
title('RGB (red-green-blue) painted')
hp1=gpatch(F,V,C_rgb_1);
hp2=gpatch(F,V2,C_rgb_2);
hp3=gpatch(F,V3,C_rgb_3);
legend([hp1 hp2 hp3],{'RGB colored with map 1','RGB colored with map 2','RGB colored with map 3'});
axisGeom;
camlight headlight; 

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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
