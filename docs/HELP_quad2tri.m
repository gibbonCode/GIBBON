%% quad2tri
% Below is a demonstration of the features of the |quad2tri| function

%% Syntax
% |[varargout]=quad2tri(varargin);|

%% Description 
% 
%% Examples 
% 

%%
clear; close all; clc;

% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=1;

%% Examples 
% 
%% Example: Converting a quandrangulated surface to a triangulated surface

%%

[X,Y,Z]=peaks(25);
[F,V]=surf2patch(X,Y,Z);


%%
% Convert triangular faces to quadrilateral faces
[Ft,Vt]=quad2tri(F,V);

%%
% Visualisation

hf=figuremax(fig_color,fig_colordef); 
subplot(1,2,1);
title('Quadrangulation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','k');
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;

subplot(1,2,2);
title('Triangulation','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',Ft,'Vertices',Vt,'FaceColor','b','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','r');
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;

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
% ********** _license boilerplate_ **********
% 
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
