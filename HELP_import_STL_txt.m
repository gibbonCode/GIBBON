%% import_STL_txt
% Below is a demonstration of the features of the |import_STL_txt| function

%%
clear; close all; clc; 

%%
% Plot settings
fig_color='w'; 
fig_colordef='white'; 
faceAlpha=1;
fontSize=10; 

%% Import STL file as patch data

loadNameIndentor='C:\Users\kmmoerman\Dropbox\DAVID_KEVIN\FEA_2014_06_18\FIT_Socket_Data_files\FITSocket_indentorstl_new.STL';

[stlStruct] = import_STL_txt(loadNameIndentor);

F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};

%%
% Plotting the models 

figuremax(fig_color,fig_colordef);
title('Imported patch data from multi-solid STL','fontSize',fontSize);
xlabel('X','fontSize',fontSize);ylabel('Y','fontSize',fontSize); zlabel('Z','fontSize',fontSize); hold on;

patch('Faces',F,'Vertices',V,'FaceColor','r','EdgeColor','k','FaceAlpha',faceAlpha);

view(3); axis equal; axis tight; axis vis3d; grid on; 
camlight('headlight');
lighting flat;
drawnow;

