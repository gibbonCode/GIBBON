%% import_STL_bin
% Below is a demonstration of the features of the |import_STL_bin| function

%% Syntax
% |[stlStruct] = import_STL_bin(fileName);|

%% Description
% Use |import_STL_bin| to import binary type STL files.

%% Examples

clear; close all; clc;

%%
% Plot settings
fontSize=25; 

%% Import a binary type STL file as patch data

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'femur_bin.stl'); 
[stlStruct] = import_STL_bin(fileName);

F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};

%%
% Merging nodes example
[~,ind1,ind2]=unique(pround(V,5),'rows');
V=V(ind1,:);
F=ind2(F);

%%
% Plotting the model

cFigure;
title('Imported patch data from STL','fontSize',fontSize);
gpatch(F,V,'g');
axisGeom;
camlight('headlight');
lighting phong; axis off; 
drawnow;


%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
