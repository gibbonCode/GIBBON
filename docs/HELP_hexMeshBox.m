%% hexMeshBox
% Below is a demonstration of the features of the |hexMeshBox| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=20;
faceAlpha1=1;
faceAlpha2=0.5;

%% CREATING A MESHED BOX
boxDim=[5 6 7];
boxEl=[4 5 6];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E=meshStruct.E;
V=meshStruct.V;
F=meshStruct.F;
Fb=meshStruct.Fb;
faceBoundaryMarker=meshStruct.faceBoundaryMarker;

%%
% Plotting model
cFigure;
subplot(1,2,1);
title('Box boundaries faces','FontSize',fontSize);
hold on;

gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha1);
% [hp]=patchNormPlot(Fb,V,1); %Display face normals

axisGeom(gca,fontSize); 
colormap(gjet(6)); icolorbar; 


subplot(1,2,2);
title('Box hexahedral mesh','FontSize',fontSize);
hold on;

gpatch(F,V,0.5*ones(1,3),'k',faceAlpha2);
patchNormPlot(F,V);
axisGeom(gca,fontSize); 
camlight headlight;

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
