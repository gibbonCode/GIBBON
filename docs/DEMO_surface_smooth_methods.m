%% DEMO_surface_smooth_methods
% Below is a demonstration for:
% 
% * The use of |patchSmooth| to smooth surface using either the Laplacian or Humphreys classes method

%%

clear; close all; clc; 

%%
fontSize=15;
cMap=gjet(250);

%%
C=gjet(4);
c1=C(4,:);
c2=C(2,:);
[F,V]=stanford_bunny; %Some graphics data

V_noisy=V+5*randn(size(V));

V1=V_noisy;
V2=V_noisy;

nSteps=15;
nSub=1;

%% Comparing smoothening methods

hf=cFigure; 

subplot(1,2,1); hold on; 
ht1=title('Laplacian smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp1=gpatch(F,V1,c1,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);

subplot(1,2,2); hold on; 
ht2=title('Humphreys classes smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp2=gpatch(F,V2,c2,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);

drawnow;

%%


animStruct.Time=linspace(0,1,nSteps);
animStruct.Handles{1}=[hp1,hp2,ht1,ht2]; %Handles of objects to animate
animStruct.Props{1}={'Vertices','Vertices','String','String'}; %Properties of objects to animate
animStruct.Set{1}={V1,V2,'Laplacian smoothing i=0','Humphreys classes smoothing i=0'}; %Property values for to set in order to animate
    
clear cPar;

cPar1.n=nSub; %Number of iterations
cPar1.Method='LAP'; %Smooth method

cPar2.n=nSub; %Number of iterations
cPar2.Method='HC'; %Smooth method

for q=2:1:nSteps    
    
    [V1]=patchSmooth(F,V1,[],cPar1);
    
    [V2]=patchSmooth(F,V2,[],cPar2);
    
    %Set entries in animation structure
    animStruct.Handles{q}=[hp1,hp2,ht1,ht2]; %Handles of objects to animate
    animStruct.Props{q}={'Vertices','Vertices','String','String'}; %Properties of objects to animate
    animStruct.Set{q}={V1,V2,['Laplacian smoothing i=',num2str(q*nSub)],['Humphreys classes smoothing i=',num2str(q*nSub)]}; %Property values for to set in order to animate
end

anim8(hf,animStruct);

hf=cFigure; 

subplot(1,2,1); hold on; 
ht1=title('Laplacian smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp1=gpatch(F,V1,c1,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);

subplot(1,2,2); hold on; 
ht2=title('Humphreys classes smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp2=gpatch(F,V2,c2,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);

drawnow;

%% Comparing using color

clear hf animStruct;
V1=V_noisy;
V2=V_noisy;

C1=minDist(V1,V);
C2=minDist(V2,V);

hf=cFigure; 

subplot(1,2,1); hold on; 
ht1=title('Laplacian smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp1=gpatch(F,V1,C1,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);
colormap(cMap); caxis manual; caxis([0 max(C1(:))]); hc=caxis;
colorbar;

subplot(1,2,2); hold on; 
ht2=title('Humphreys classes smoothing','FontSize',fontSize);
gpatch(F,V,0.5*ones(1,3),'none',0.2);
hp2=gpatch(F,V2,C2,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);
colormap(cMap); caxis(hc);
colorbar;

drawnow;

%%

animStruct.Time=linspace(0,1,nSteps);
animStruct.Handles{1}=[hp1,hp2,ht1,ht2,hp1,hp2]; %Handles of objects to animate
animStruct.Props{1}={'Vertices','Vertices','String','String','CData','CData'}; %Properties of objects to animate
animStruct.Set{1}={V1,V2,'Laplacian smoothing i=0','Humphreys classes smoothing i=0',C1,C2}; %Property values for to set in order to animate
    
clear cPar;

cPar1.n=nSub; %Number of iterations
cPar1.Method='LAP'; %Smooth method

cPar2.n=nSub; %Number of iterations
cPar2.Method='HC'; %Smooth method

for q=2:1:nSteps    
    
    [V1]=patchSmooth(F,V1,[],cPar1);    
    [V2]=patchSmooth(F,V2,[],cPar2);
    
    C1=minDist(V1,V);
    C2=minDist(V2,V);
    
    %Set entries in animation structure
    animStruct.Handles{q}=[hp1,hp2,ht1,ht2,hp1,hp2]; %Handles of objects to animate
    animStruct.Props{q}={'Vertices','Vertices','String','String','CData','CData'}; %Properties of objects to animate
    animStruct.Set{q}={V1,V2,['Laplacian smoothing i=',num2str(q*nSub)],['Humphreys classes smoothing i=',num2str(q*nSub)],C1,C2}; %Property values for to set in order to animate
end

anim8(hf,animStruct);

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
