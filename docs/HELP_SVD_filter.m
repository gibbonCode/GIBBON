%% SVD_filter
% Below is a demonstration of the features of the |SVD_filter| function

%%
clear; close all; clc;

%% Syntax
% |[Zm]=SVD_filter(Z,P,T);|

%% Description 
% This funciton uses singular value decomposition to smoothen 2D matrix
% data

%% Examples 
% 

%%
% Plot settings
font_size=20;
cmap=resampleColormap([1 0 0; 1 1 1; 0 0 1],250);

%%
% Create example data 

%Create clean data
n=35;
s=5;
[X,Y]=ndgrid(linspace(-4*s,4*s,n));
Z=n*exp( -0.5.*((X./s).^2+(Y./s).^2));
Z(X<0)=-Z(X<0); %Add sharp feature

%Create noise eroded data
Zn=Z+n/30*randn(size(Z)); 

%% Using |SVD_filter|
P=[1-1e-4 1e-4]; 
T=1;
[Zm]=SVD_filter(Zn,P,T);

[F,V]=surf2patch(Z);
[~,Vn]=surf2patch(Zn);
[~,Vm]=surf2patch(Zm);
Cn=Vn(:,3)-V(:,3);
Cm=Vm(:,3)-V(:,3);

c=max(abs(Cn(:)));

%%
% Visualize

cFigure;
subplot(1,3,1); hold on;
title('Clean','FontSize',font_size);
hp=gpatch(F,V,V(:,3),'k');
hp.FaceColor='interp';
axisGeom(gca,font_size);
colormap(gca,gjet(250));
camlight headlight;

subplot(1,3,2); hold on;
title('Raw','FontSize',font_size);
hp=gpatch(F,Vn,Cn,'k');
hp.FaceColor='interp';
axisGeom(gca,font_size);
colormap(gca,cmap); colorbar;
caxis([-c,c]);
camlight headlight;

subplot(1,3,3); hold on
title('SVD filtered','FontSize',font_size);
% gpatch(F,V,'g','none',0.5);
hp=gpatch(F,Vm,Cm,'k');
hp.FaceColor='interp';
axisGeom(gca,font_size);
colormap(gca,cmap); colorbar;
camlight headlight;
caxis([-c,c]);
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
