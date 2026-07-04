%% TriScatteredInterp_ND
% Below is a demonstration of the features of the |TriScatteredInterp_ND| function

%%
clear; close all; clc;

%% Syntax
% |[DI]=TriScatteredInterp_ND(DT,D,XI,InterpMethod);|
% |[DI]=TriScatteredInterp_ND(P,D,XI,InterpMethod);|

%% Description 
% DEPRICATED
% See also: |scatteredInterpolant|

%% Examples 
% 

%%
% Plot settings
cMap=viridis(250);
fontSize=15;

%%

[X,Y,Z]=ndgrid(linspace(0,2*pi,10)); %Perfect grid
P=[X(:) Y(:) Z(:)]; %Grid as nx3 array
P=P+0.1*randn(size(P)); %Add noise so it is not a perfect grid

%Compute multidimensional (here 3D vector) quantity
UX=cos(P(:,3));
UY=0.25*P(:,1)./(2*pi);
UZ=cos(P(:,1));

U=[UX UY UZ];
C=sqrt(sum(U.^2,2)); %Vector lenghts

%%

[Xi,Yi,Zi]=ndgrid(linspace(0,2*pi,20)); % Interpolation grid
Pi=[Xi(:) Yi(:) Zi(:)]; %Interpolation set as nx3 array

DT=delaunayTriangulation(P);
InterpMethod='natural';
Ui=TriScatteredInterp_ND(DT,U,Pi,InterpMethod);
Ci=sqrt(sum(Ui.^2,2)); %Vector lenghts
%%

cFigure;
subplot(1,2,1); hold on; 
title('Input');
quiverVec(P,U,0.5,C);
colormap(cMap); colorbar;
axisGeom(gca,fontSize);
camlight headlight; lighting flat

subplot(1,2,2); hold on; 
title('Interpolated');
quiverVec(Pi,Ui,0.5,Ci);
colormap(cMap); colorbar;
axisGeom(gca,fontSize);
camlight headlight; lighting flat

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
