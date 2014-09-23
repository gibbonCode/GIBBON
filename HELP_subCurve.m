%% subCurve
% Below is a demonstration of the features of the |subCurve| function
%%

clear; close all; clc;

%%
% PLOT SETTINGS
figColor='w'; figColorDef='white';
fontSize=15;
markerSize1=45;
lineWidth1=2;
lineWidth2=5;
lineWidth3=2;
faceAlpha=0.5;

%% LINEARLY UPSAMPLING A CURVE

%%
% Simulating a curve
Vt=[0 0 0; 10 0 0; 5 10 0; 10 0 10; 0 10 10; ];

%Setting number of desired intermediate points to be added 
np=3; 

[VN]=subCurve(Vt,3);

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);
title('A linearly upsampled curve','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Vt,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(VN,'r.-','lineWidth',lineWidth1/2,'MarkerSize',markerSize1/2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)

%%