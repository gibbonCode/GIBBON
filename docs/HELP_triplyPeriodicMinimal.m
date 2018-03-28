%% triplyPeriodicMinimal
% Below is a demonstration of the features of the |triplyPeriodicMinimal| function

%%
clear; close all; clc;

%%
% Plot settings
cMap=jet(250);
faceAlpha1=1;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';
fontSize=25; 

%% SURFACE VISUALIZATIONS

pColors=gjet(6);

n=50;

[X,Y,Z]=meshgrid(linspace(-pi,pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'p');
[F,V] = isosurface(X,Y,Z,S,0.1);

%%
% Visualization

cFigure;

subplot(2,3,1);
title('Schwarz P-surface','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(1,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

[X,Y,Z]=meshgrid(linspace(-pi,pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'d');
[F,V] = isosurface(X,Y,Z,S,0.1);

subplot(2,3,2);
title('d','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(2,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

[X,Y,Z]=meshgrid(linspace(-2*pi,2*pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'g');
[F,V] = isosurface(X,Y,Z,S,0.6);

subplot(2,3,3);
title('Gyroid','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(3,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

[X,Y,Z]=meshgrid(linspace(-pi,pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'n');
[F,V] = isosurface(X,Y,Z,S,0);

subplot(2,3,4);
title('Neovius','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(4,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

[X,Y,Z]=meshgrid(linspace(-pi,pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'w');
[F,V] = isosurface(X,Y,Z,S,-0.1);

subplot(2,3,5);
title('w','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(5,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

[X,Y,Z]=meshgrid(linspace(-pi,pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'pw');
[F,V] = isosurface(X,Y,Z,S,0.5);

subplot(2,3,6);
title('pw','FontSize',fontSize);
hold on;
gpatch(F,V,pColors(6,:),'none',1);
axisGeom; 
camlight headlight;
view(-50,30);

drawnow;

%% EXAMPLE SIMULATING REGULARIZED TRABECULAE OR POROUS MEDIA

n=50;
[X,Y,Z]=meshgrid(linspace(-2*pi,2*pi,n));
S=triplyPeriodicMinimal(X,Y,Z,'g');

T=0.6;      
c=(max(S(:))*T);
L=S>=c;
[F1,V1,C1]=im2patch(S,L,'vb'); 

L=S<c;
[F2,V2,C2]=im2patch(S,L,'vb'); 

%%
% Visualization

cFigure;
subplot(1,2,1);
title('Gyroid based "network"','FontSize',fontSize);
xlabel('J','FontSize',fontSize);ylabel('I','FontSize',fontSize); zlabel('K','FontSize',fontSize); hold on;
gpatch(F1,V1,C1,'k',1);
axisGeom; 
colormap(cMap); 
caxis([min(S(:)) max(S(:))]);
camlight headlight;

subplot(1,2,2);
title('Gyroid based "matrix"','FontSize',fontSize);
xlabel('J','FontSize',fontSize);ylabel('I','FontSize',fontSize); zlabel('K','FontSize',fontSize); hold on;
gpatch(F2,V2,C2,'k',1);
axisGeom; 
colormap(cMap); 
caxis([min(S(:)) max(S(:))]);
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
