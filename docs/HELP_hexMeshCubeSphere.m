%% hexMeshCubeSphere
% Below is a demonstration of the features of the |hexMeshCubeSphere| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=15;
faceAlpha1=0.5;

%%

cPar.boxWidth=1;
cPar.OuterSphereRadius=cPar.boxWidth/6; 
cPar.InnerSphereRadius=cPar.OuterSphereRadius/2; 
cPar.CoreSphereRadius=cPar.InnerSphereRadius/2; 
cPar.numElementsCube=16;
cPar.numElementsCubeSphere=6;
cPar.numElementsSphereMantel=6;
cPar.numElementsSphereCore=6;
cPar.nSmooth=15;

[E,V,CE,Fb,Cb]=hexMeshCubeSphere(cPar);

%%

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L,:),'hex8');

%%

cFigure;
subplot(1,2,1); hold on;
title('Cut-view of the mesh','FontSize',fontSize);

gpatch(Fs,V,Cs);

axisGeom(gca,fontSize);
colormap(gca,gjet(3)); icolorbar; 
camlight headlight;

subplot(1,2,2); hold on;
title('Mesh boundaries','FontSize',fontSize);

gpatch(Fb,V,Cb,'none',0.5);

axisGeom(gca,fontSize);
colormap(gca,gjet(8)); icolorbar; 
camlight headlight;

drawnow; 

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
