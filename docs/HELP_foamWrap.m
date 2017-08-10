%% foamWrap
% Below is a demonstration of the features of the |foamWrap| function

%% Syntax
% |[FT,VT,CT,CT_c]=foamWrap(F,V,C,cPar);|

%% Description
% Use |foamWrap| to generate a foam like structure on top of an input mesh

%% Examples

clear; close all; clc;

%% 
% Plot Settings

fontSize=15;
faceAlpha=1;
edgeColor=0.1*ones(1,3);
edgeWidth=1;
cmap=gjet(250);

%% 
% Create surface model 

[F,V,~]=geoSphere(2,1); %Geodesic sphere
% [F,V]=parasaurolophus;
% [F,V]=cow;
% [F,V]=stanford_bunny;

[F,V,C,indIni]=triPolyDualRefine(F,V);

%%
cFigure;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','k');

camlight headlight; 
colormap(cmap);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on; axis off; 
view(0,59);

%%

cPar.n=3; 
cPar.dirFlip=1; 
cPar.foamThickness=[]; %Empty uses default which is mean edgelength based
cParSmooth.Method='HC';
cParSmooth.n=25;
cPar.cParSmooth=cParSmooth; 

%%
L_remove=true(size(F,1),1);
[FT,VT,CT,CT_c]=foamWrap(F,V,C,cPar);

%%
cFigure; hold on; 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

gpatch(FT,VT,CT_c,'none',1);

colormap(gray(250)); icolorbar;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
view(0,58.25);
camlight headlight; 
axis off;

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
