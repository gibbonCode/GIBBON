%% trisurf_intersect
% Below is a demonstration of the features of the |trisurf_intersect| function

%%
clear; close all; clc;

%% Syntax
% |[P]=trisurf_intersect(TRI1c,N1c,Vn1c,V1,TRI2c,N2c,Vn2c,V2,TRI3c,N3c,Vn3c,V3,d);|

%% Description 
% This function computes the intersection of a triplet of triangulated
% surfaces. 

%% Examples 
% 

%Surface 1
n1=10;
[X1,Y1]=meshgrid(linspace(-3,3,n1));
Z1=peaks(X1,Y1); 
Z1=Z1./max(abs(Z1(:)));
[F1,V1]=grid2patch(X1,Y1,Z1);
[F1,V1]=quad2tri(F1,V1,'a');
[N1,Vn1]=patchNormal(F1,V1);

n2=10;
[X2,Y2]=meshgrid(linspace(-3,3,n2));
Z2=peaks(X2,Y2); 
Z2=Z2./max(abs(Z2(:)));
[F2,V2]=grid2patch(X2,Y2,Z2);
[F2,V2]=quad2tri(F2,V2,'a');
R2=euler2DCM([0.4*pi,0,0]);
V2=V2*R2;
[N2,Vn2]=patchNormal(F2,V2);

n3=10;
[X3,Y3]=meshgrid(linspace(-3,3,n3));
Z3=peaks(X3,Y3); 
Z3=Z3./max(abs(Z3(:)));
[F3,V3]=grid2patch(X3,Y3,Z3);
[F3,V3]=quad2tri(F3,V3,'a');
R3=euler2DCM([0,0.45*pi,0]);
V3=V3*R3;
[N3,Vn3]=patchNormal(F3,V3);

%%

[P]=trisurf_intersect(F1,N1,Vn1,V1,F2,N2,Vn2,V2,F3,N3,Vn3,V3,5);

%%

cFigure; hold on; 
gpatch(F1,V1,'rw','k',0.5);
gpatch(F2,V2,'gw','k',0.5);
gpatch(F3,V3,'bw','k',0.5);
plotV(P,'k.','markerSize',50)
axisGeom; camlight headlight; 
gdrawnow

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
