%% curvePathOrderFix
% Below is a demonstration of the features of the |curvePathOrderFix| function

%%
clear; close all; clc;

%% Syntax
% |[Vf,indFix]=curvePathOrderFix(V);|

%% Description 
% This function unscrambles a scrambles curve using distances. The function
% assumes that the proper point order can be resolved by choosing the next
% nearest point. 

%% Examples 
% 

t=linspace(0,2*pi,25)';
V=[cos(t) sin(t) zeros(size(t))];
V=V(randperm(size(V,1)),:); %Scramble point set

%% 
% Undo scrambling using |curvePathOrderFix|

[Vf,indFix]=curvePathOrderFix(V);

%%
% Visualization

cFigure; 
subplot(1,2,1); hold on;
title('Scrambled curve');
plotV(V,'r.-','LineWidth',2,'MarkerSize',25);
axis tight; axis equal; grid on; box on; view(2); 

subplot(1,2,2); hold on;
title('Unscrambled curve');
plotV(Vf,'g.-','LineWidth',2,'MarkerSize',25);
axis tight; axis equal; grid on; box on; view(2); 
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
