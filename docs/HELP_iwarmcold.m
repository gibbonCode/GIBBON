%% iwarmcold
% Below is a demonstration of the features of the |iwarmcold| function

%%
clear; close all; clc;

%% Syntax
% |[cMap]=iwarmcold(n);|

%% Description 
% The iwarmcold colormap

%% Examples 
% 

%%
%Plot settings
fontSize=15;

%% Example 1: Accessing the colormap

cMap=iwarmcold(250); %Get 25 color levels from the colormap

%% Example 2: Applying/using the colormap

% Create example data for visualizations
n=250;
s=1;
[X,Y]=ndgrid(linspace(-3*s,3*s,n));
Z=exp( -0.5.*((X./s).^2+(Y./s).^2));
Z=Z./max(Z(:));
Z(X<0)=-Z(X<0);
colorLim=[-1 1];

%%
cFigure; hold on;
title('gjet colormap','FontSize',fontSize);
imagesc(Z);
colormap(cMap); colorbar;
axis tight; axis equal; axis xy; box on;
axis off;
set(gca,'FontSize',fontSize);
clim(colorLim);
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
