%% import_off
% Below is a demonstration of the features of the |import_off| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; 
fig_colordef='white'; 
faceAlpha=1;
fontSize=10; 

%% Import OFF file as patch data

%Set main folder
defaultFolder = fileparts(mfilename('fullpath'));
pathName=fullfile(defaultFolder,'data','OFF'); 
testCase=1; 
switch testCase
    case 1
        offName='elephant-50kv.off';
end
fileName=fullfile(pathName,offName); 

%%
% The |import_off| function imports <http://www.geomview.org> type .off
% files specifying a surface mesh (provided all faces are of the same
% type!). The output consists of (patch data) faces F and vertices V. 

[F,V] = import_off(fileName);

%%
% Plotting the models 

figuremax(fig_color,fig_colordef);
title('Imported patch data from OFF file','fontSize',fontSize);
xlabel('X','fontSize',fontSize);ylabel('Y','fontSize',fontSize); zlabel('Z','fontSize',fontSize); hold on;

patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'EdgeColor','k','FaceAlpha',faceAlpha);

view(3); axis equal; axis tight; axis vis3d; grid on; axis off; 
view(194.5,-40);
camlight('headlight');
lighting flat;
drawnow;

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
