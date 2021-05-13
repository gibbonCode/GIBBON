%% diamondLattice
% Below is a demonstration of the features of the |diamondLattice| function

%%
clear; close all; clc;

%% Syntax
% |[Ep,Et,VT,Ct]=diamondLattice(sampleSize,nRepeat,strutThickness,plotOn);|

%% Description 
% UNDOCUMENTED 

%%
% Plotting settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=25;
markerSize2=10;
cMap=gjet(4);

%% Examples 
% 

%Latticeparameters
nRepeat=3; %Number of repetitions of the lattice pattern
sampleSize=30;
nSubPenta=2;
strutThickness=1; %Set the strut thickness

%% Create diamond lattice

[Ep,Et,VT,Ct]=diamondLattice(sampleSize,nRepeat,strutThickness,0);
[Ep,VT]=subPenta(Ep,VT,nSubPenta,3); %Sub-divide pentahedra
% strutThicknessCheck=mean(patchEdgeLengths(Fp{1},VT));

%Get element faces for visualization
Fp=element2patch(Ep,[],'penta6');
Ft=element2patch(Et,[],'tet4');

%% Visualization

cFigure; 
subplot(1,2,1); hold on; 
gpatch(Fp,VT,'bw','none',1);
gpatch(Ft,VT,'bw','none',1);
axisGeom; camlight headlight; 

subplot(1,2,2); hold on; 
hpl=gpatch(Fp,VT,'rw','r',0.5);
hpl(end+1)=gpatch(Ft,VT,'gw','g',0.5);
legend(hpl,{'Pentahedral triangles','Pentahedra quads','Tetrahedral triangles'});
axisGeom; camlight headlight; 

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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
