%% rhombicDodecahedronMesh
% Below is a demonstration of the features of the |rhombicDodecahedronMesh| function

%%
clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceAlpha1=1;

%% Creating a mesh of rhombic dodecahedra

r=sqrt(2)/2; %Radii, the chosen level results in X,Y spacings of 1
n=3; %Desired number of copies in each direction 

%The actual input 
nCopies=[n+1 n+1 n+ceil((n+1)/2)]; %Number of offset copies

[Fc_Q,Fc_T,Ft_Q,Ft_T,Ct_Q,Ct_T,Vt]=rhombicDodecahedronMesh(r,nCopies);

%%
% Plotting results
cFigure; hold on;
title('A mesh of rhombicDodecahedra','FontSize',fontSize);
gpatch(Ft_Q,Vt,Ct_Q,'k',faceAlpha1);
colormap(gjet);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow; 

%%

cFigure; hold on;
title('A mesh of rhombicDodecahedra plotted with triangular faces','FontSize',fontSize);
gpatch(Ft_T,Vt,Ct_T,'k',faceAlpha1);
colormap(gjet);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow; 

%% EXAMPLE: Creating a "bubble image"
% The triangular surface output can but used with the function |triSurf2Im|
% to convert the patch data to image data where the interior and exterior
% are labelled differently. For the mesh created the result is like a
% highly regular bubble image. 
voxelSize=0.05;
[M,G,ML]=triSurf2Im(Ft_T,Vt,voxelSize);
imOrigin=G.origin;

%%
% Plot results

cFigure;
hold on;
% title('"Bubble image"','FontSize',fontSize);

gpatch(Ft_Q,Vt,0.5*ones(1,3),'none',0.1);

L_plot=false(size(ML));
L_plot(:,:,round(size(ML,3)/2))=1;
L_plot=L_plot&~isnan(ML);
[Fm,Vm,Cm]=ind2patch(L_plot,ML,'sk');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
gpatch(Fm,Vm,Cm,'none',faceAlpha1);

L_plot=false(size(ML));
L_plot(round(size(ML,1)/2),:,:)=1;
L_plot=L_plot&~isnan(ML);
[Fm,Vm,Cm]=ind2patch(L_plot,ML,'si');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
gpatch(Fm,Vm,Cm,'none',faceAlpha1);

L_plot=false(size(ML));
L_plot(:,round(size(ML,2)/2),:)=1;
L_plot=L_plot&~isnan(ML);
[Fm,Vm,Cm]=ind2patch(L_plot,ML,'sj');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
gpatch(Fm,Vm,Cm,'none',faceAlpha1);

colormap(gjet(max(ML(:)))); colorbar;
camlight('headlight'); 
axisGeom(gca,fontSize);

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
