%% import_STL_txt
% Below is a demonstration of the features of the |import_STL_txt| function

%%
clear; close all; clc; 

%%
% Plot settings
fig_color='w'; 
fig_colordef='white'; 
faceAlpha=0.5;
fontSize=10; 

%% Import STL file as patch data

%Set main folder
defaultFolder = fileparts(mfilename('fullpath'));
pathName=fullfile(defaultFolder,'data','STL'); 
testCase=1; 
switch testCase
    case 1
        stlName='standford_bunny_multi.stl';
    case 2
        stlName='femur.stl';
    case 3
        stlName='hip_implant.stl';
end
fileName=fullfile(pathName,stlName); 
[stlStruct] = import_STL_txt(fileName);

%%
% Plotting the models 
pColors=autumn(numel(stlStruct.solidNames));

figuremax(fig_color,fig_colordef);
title('Imported patch data from multi-solid STL','fontSize',fontSize);
xlabel('X','fontSize',fontSize);ylabel('Y','fontSize',fontSize); zlabel('Z','fontSize',fontSize); hold on;
for q=1:1:numel(stlStruct.solidNames)
    F=stlStruct.solidFaces{q};
    V=stlStruct.solidVertices{q};
    patch('Faces',F,'Vertices',V,'FaceColor',pColors(q,:),'EdgeColor','k','FaceAlpha',faceAlpha);
end
view(3); axis equal; axis tight; axis vis3d; grid on; 
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