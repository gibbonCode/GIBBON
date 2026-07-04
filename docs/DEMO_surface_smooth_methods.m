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

[F,V]=stanford_bunny; %Some graphics data

V_noisy=V+5*randn(size(V));

V1=V_noisy;
V2=V_noisy;

nSteps=25;

%% Comparing smoothening methods

hf=cFigure; 

subplot(1,2,1); hold on; 
ht1=title('Laplacian smoothing','FontSize',fontSize);
gpatch(F,V,'kw','none',0.2);
hp1=gpatch(F,V1,'rw','r',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
view(0,0);

subplot(1,2,2); hold on; 
ht2=title('Humphreys classes smoothing','FontSize',fontSize);
gpatch(F,V,'kw','none',0.2);
hp2=gpatch(F,V2,'gw','g',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
view(0,0);

drawnow;

%%

stepRange=0:1:nSteps;  
animStruct.Time=stepRange;
animStruct.Handles{1}=[hp1,hp2,ht1,ht2]; %Handles of objects to animate
animStruct.Props{1}={'Vertices','Vertices','String','String'}; %Properties of objects to animate
animStruct.Set{1}={V1,V2,'Laplacian smoothing i=0','Humphreys classes smoothing i=0'}; %Property values for to set in order to animate
    
cPar1.Method='LAP'; %Smooth method
cPar2.Method='HC'; %Smooth method

q=1;
for c=stepRange      
    cPar1.n=c; %Number of iterations
    [V1s]=patchSmooth(F,V1,[],cPar1);
    
    cPar2.n=c; %Number of iterations
    [V2s]=patchSmooth(F,V2,[],cPar2);
    
    %Set entries in animation structure
    animStruct.Handles{q}=[hp1,hp2,ht1,ht2]; %Handles of objects to animate
    animStruct.Props{q}={'Vertices','Vertices','String','String'}; %Properties of objects to animate
    animStruct.Set{q}={V1s,V2s,['Laplacian smoothing i=',num2str(c)],['Humphreys classes smoothing i=',num2str(c)]}; %Property values for to set in order to animate
    q=q+1;
end

anim8(hf,animStruct);

%% Compute distance metric to original

C1=minDist(V1s,V);
C2=minDist(V2s,V);

hf=cFigure; 

subplot(1,2,1); hold on; 
ht1=title('Laplacian smoothing','FontSize',fontSize);
gpatch(F,V,'kw','none',0.2);
hp1=gpatch(F,V1s,C1,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);
colormap(cMap); caxis manual; caxis([0 max(C1(:))]); hc=caxis;
colorbar;

subplot(1,2,2); hold on; 
ht2=title('Humphreys classes smoothing','FontSize',fontSize);
gpatch(F,V,'kw','none',0.2);
hp2=gpatch(F,V2s,C2,'k',1);

axisGeom(gca,fontSize); axis manual; axis off; 
camlight right;
zoom(1.5);
colormap(cMap); caxis(hc);
colorbar;

drawnow

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
