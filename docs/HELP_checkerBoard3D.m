%% checkerBoard3D
% Below is a demonstration of the features of the |checkerBoard3D| function

%%
clear; close all; clc;

%% Syntax
% |M=checkerBoard3D(siz);|

%% Description 
% This function creates a checkboard image of the size siz whereby elements
% are either black (0) or white (1). The first element is white.

%% Examples 

%%
% Plot settings
fontSize=15; 

%% Example creating a 2D checkerboard image

siz=[4 6]; %Image size
M=checkerBoard3D(siz); %Create checkerboard image

%%
% Plotting results

cFigure;
title('A 2D checkerboard pattern','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
imagesc(M);

colormap gray;
axis equal; view(2); axis tight;
set(gca,'FontSize',fontSize);
drawnow;

%% Example creating a 3D checkerboard image

siz=[6 6 3]; %Image size
M=checkerBoard3D(siz); %Create checkerboard image

%%
% Plotting results

[Fv,Vv,Cv]=ind2patch(1:numel(M),M,'vu'); %Create patch data for plotting

cFigure;
title('A 3D checkerboard pattern','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hp= patch('Faces',Fv,'Vertices',Vv,'FaceColor','flat','CData',Cv,'EdgeColor','r','FaceAlpha',1);

camlight('headlight');
colormap gray;
axis equal; view(3); axis tight;  grid on; box on; 
set(gca,'FontSize',fontSize);
drawnow;

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
