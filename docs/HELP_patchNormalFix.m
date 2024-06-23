%% patchNormalFix
% Below is a demonstration of the |patchNormalFix| function

%% Syntax
% |[F]=patchNormalFix(F);|

%% Description
% This function reorients all faces coherently such that the face normals
% point in a coherent fashion.

%% Examples

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=1;

%% Example 1: Fix face normals using |patchNormalFix|

%%
% Create example geometry with incoherent face normals

%An example surface
% [F,V]=stanford_bunny('g');
[F,V]=geoSphere(2,1);

% [F,V]=tri2quad(F,V); %Conver to quadrilaterals for testing

%Alter face orientations
VF=patchCentre(F,V);
logicFlip=VF(:,1)>0 | VF(:,3)>0;
F(logicFlip,:)=fliplr(F(logicFlip,:));
F=fliplr(F);

%%
% Fix face normals using |patchNormalFix|
[F_fix,L]=patchNormalFix(F);

%%
% Visualisation
cFigure;
subplot(1,2,1); hold on; 
title('Incoherent faces','FontSize',fontSize);
gpatch(F,V,'g');
gpatch(F(1,:),V,'none','r',1,3);
patchNormPlot(F,V);
axisGeom(gca,fontSize);
camlight('headlight');

subplot(1,2,2); hold on; 
title('Coherent faces','FontSize',fontSize);
gpatch(F_fix,V,L);
patchNormPlot(F_fix,V);
axisGeom(gca,fontSize);
colormap(gjet(2)); icolorbar; 
camlight('headlight');

drawnow; 

%% Example 2: Control face orientation by providing a proper face

indStart=250; %Index of face with correct orientation
[F_fix,L]=patchNormalFix(F,indStart);

%%
% Visualisation
cFigure;
subplot(1,2,1); hold on; 
title('Incoherent faces','FontSize',fontSize);
gpatch(F,V,'g');
gpatch(F(indStart,:),V,'none','r',1,3);
patchNormPlot(F,V);
axisGeom(gca,fontSize);
camlight('headlight');

subplot(1,2,2); hold on; 
title('Coherent faces','FontSize',fontSize);
gpatch(F_fix,V,L);
patchNormPlot(F_fix,V);
axisGeom(gca,fontSize);
colormap(gjet(2)); icolorbar; 
camlight('headlight');

drawnow; 

%% Example 3: Control face orientation by providing a proper face

logicStart=~logicFlip; 
[F_fix,L]=patchNormalFix(F,logicStart);

%%
% Visualisation
cFigure;
subplot(1,2,1); hold on; 
title('Incoherent faces','FontSize',fontSize);
gpatch(F,V,'g');
gpatch(F(logicStart,:),V,'none','r',1,3);
patchNormPlot(F,V);
axisGeom(gca,fontSize);
camlight('headlight');

subplot(1,2,2); hold on; 
title('Coherent faces','FontSize',fontSize);
gpatch(F_fix,V,L);
patchNormPlot(F_fix,V);
axisGeom(gca,fontSize);
colormap(gjet(2)); icolorbar; 
camlight('headlight');

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
