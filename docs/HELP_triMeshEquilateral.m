%% triMeshEquilateral
% Below is a demonstration of the features of the |triMeshEquilateral| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
font_size=20;
cmap=gray(250);
falpha=1;
patch_types={'sx','sy','sz','v'};
ptype=3;
no_slices=4;
mark_siz1=25;
mark_siz2=25;
mark_siz3=15;
line_width1=2;
F_alpha1=1;
F_alpha2=0.8;

%% 
% Control parameters

%Desired mesh point spacing
pointSpacing=1;

%Mesh region extrema
maxV=[5 6];
minV=[-5 -5];

%% CREATING AN EQUILATERAL TRIANGLE MESH

[F,V]=triMeshEquilateral(minV,maxV,pointSpacing);

%%
% Plottting mesh
hf1=cFigure;
title('Equilateral mesh','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size); zlabel('Z','FontSize',font_size);
hold on;
hpm=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','g','FaceAlpha',F_alpha1,'lineWidth',line_width1);
% [hp]=patchNormPlot(F,V,1);
colormap autumn; 
axis equal; view(2); axis tight;  grid on; 
set(gca,'FontSize',font_size);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
