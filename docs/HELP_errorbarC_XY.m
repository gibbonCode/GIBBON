%% errorbarC_XY
% Below is a demonstration of the features of the |errorbarC_XY| function

%%
clear; close all; clc;

%% Syntax
% |H=errorbarC_XY(x,y,em,ep,w,C,d);|

%% Description 
% Visualizes colormapped arror bars. 

%% Examples 
%

%%
%
cMap=gjet(25);
markerSize=50; 
lineWidth1=3;
lineWidth2=5;

%%
n=25;
x=linspace(0,2*pi,n)';
y=sin(x);

eUpper=0.05+0.5*rand(size(x)); %Lower deviation
eLower=eUpper; %Upper deviation
w=0.2; %Width 

cData=eUpper;
C=cmaperise(cData,gjet(n),[min(cData) max(cData)]);

%%

d=1; 

%%

cFigure; hold on; 
plot(x,y,'k.-','MarkerSize',markerSize,'LineWidth',lineWidth1);

H=errorbarC_XY(x,y,eUpper,eLower,w,C,d);
set(H,'LineWidth',lineWidth2);

axis tight; axis equal; box on; grid on; 
set(gca,'FontSize',20); 
colormap(cMap)
drawnow; 

%%

d=2; 

%%

cFigure; hold on; 
plot(x,y,'k.-','MarkerSize',markerSize,'LineWidth',lineWidth1);

H=errorbarC_XY(x,y,eUpper,eLower,w,C,d);
set(H,'LineWidth',lineWidth2);

axis tight; axis equal; box on; grid on; 
set(gca,'FontSize',20); 
colormap(cMap)
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
