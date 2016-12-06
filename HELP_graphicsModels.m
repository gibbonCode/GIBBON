%% graphicsModels
% Below is a demonstration of the features of the |graphicsModels| function

%% Syntax
% |[F,V]=graphicsModels(modelId);|

%% Description 
% 
%% Examples 

%%
clear; close all; clc;

% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=0.5;

%% 

hf=figuremax(fig_color,fig_colordef); 

cMap=gjet(7);
cNames={'Stanford bunny','Utah teapot','cow','parasaurolophus','femur','hip implant','elephant'};
for q=1:1:7
   subplot(3,3,q); 
    [F,V]=graphicsModels(q);
    
    title(cNames{q},'FontSize',fontSize);
    
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    
    hp=patch('Faces',F,'Vertices',V,'FaceColor',cMap(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','none');
    
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  axis vis3d; axis off;
    camlight('headlight'); lighting flat;
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
