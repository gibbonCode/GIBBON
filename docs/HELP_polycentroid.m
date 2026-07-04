%% polycentroid
% Below is a demonstration of the features of the |polycentroid| function

%%
clear; close all; clc;

%% Syntax
% |[Xc,Yc]=polycentroid(X,Y);|

%% Description 
% This function computes the centroid of a polygon. 
%
% Assumes row vectors or matrices whereby each row describes a polygon with
% points appearing in the order defining the polygon

%% Examples 
% 

V=[0 0 0; 1 0 0; 1 1 0; 0 1 0]; 

[Xc,Yc]=polycentroid(V(:,1)',V(:,2)')

%%

F=[1 2 3 4]; %Face

cFigure; hold on; 
gpatch(F,V,'w','k',0.5);
plotV(V,'k.','MarkerSize',25)
plot(Xc,Yc,'r.','MarkerSize',50)

axisGeom; 
gdrawnow; 

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
