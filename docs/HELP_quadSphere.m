%% quadSphere
% Below is a demonstration of the features of the |quadSphere| function

%%
clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceAlpha=0.75;
edgeColor=0.3*ones(1,3);
edgeWidth=1.5;

%% Building a quadrangulated sphere
% The function inputs are n and r which define the mesh refinement and
% radius respectively. The mesh refinement number n defines the number of
% subquadrangulation (see function |subQuad|) iterations performed on an
% initial solid. The function outputs the faces (F) and vertices (V). The
% default initial solid is the cube. However this can be altered by varying
% the optional 3rd input oriOpt. Its value sets the initial solid as:
% oriOpt=1 -> Tetrahedron
% oriOpt=2 -> Cube
% oriOpt=3 -> Octahedron
% oriOpt=4 -> Icosahedron
% oriOpt=5 -> Rhombic dodecahedron
% oriOpt=6 -> (Buckminster-Fuller) geoSphere
% For option 6, geoSphere the subdevisions are based on subtriangulations
% of the icosahedron (see |geoSphere|) after which all faces are converted
% to quadrilateral faces (see |tri2quad|).
% Below is a visualisation for n=0:1:3 for the cube as initial solid. 

% Open figure for plotting
hf=cFigure; 

%Defining geodesic dome
r=1; %sphere radius
n=0:1:3; %Refinements   
pColors=autumn(numel(n));
for q=1:1:numel(n);
    [F,V]=quadSphere(n(q),r,2);
   
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' refinement iterations'],'FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hp=patch('Faces',F,'Vertices',V);
    set(hp,'FaceColor',pColors(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    camlight headlight; 
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  grid on;
end

%% Demo showing initial shapes and subquadrangulated results

hf=cFigure; 

[F,V]=quadSphere(0,r,2);   
subplot(2,3,1); hold on;
title(['Cube 0 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

[F,V]=quadSphere(3,r,2);   
subplot(2,3,4); hold on;
title(['Cube 3 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

[F,V]=quadSphere(0,r,5);   
subplot(2,3,2); hold on;
title(['Rhombic dodecahedron 0 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','b','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

[F,V]=quadSphere(3,r,5);   
subplot(2,3,5); hold on;
title(['Rhombic dodecahedron 3 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','b','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

[F,V]=quadSphere(0,r,6);   
subplot(2,3,3); hold on;
title(['geoSphere 0 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

[F,V]=quadSphere(2,r,6);   
subplot(2,3,6); hold on;
title(['geoSphere 2 ref. it.'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

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
