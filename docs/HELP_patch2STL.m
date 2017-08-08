%% patch2STL
% Below is a demonstration of the features of the |patch2STL| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white'; 
edgeColor1='none';
edgeColor2='none';

%% 
% Get patch data
[F,V]=stanford_bunny;

%%
% Plotting the model 

figuremax(fig_color,fig_colordef);
xlabel('X');ylabel('Y'); zlabel('Z'); hold on;

patch('Faces',F,'Vertices',V,'FaceColor','b','EdgeColor','k','FaceAlpha',1);

view(3); axis equal; axis tight; axis vis3d; grid on; 
camlight('headlight');
lighting phong; axis off; 
drawnow;

%% Exporting an STL file from patch data

%Set main folder and fileName
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'standford_bunny.stl'); 

patch2STL(fileName,V,F,[],'standford_bunny');

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
