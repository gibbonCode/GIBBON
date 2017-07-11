%% hex2tet
% Below is a demonstration of the features of the |hex2tet| function

%% Syntax
% |[TET,Vtet,C]=hex2tet(HEX,V,C,tetOpt);|

%% Description 
% 

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15; 
faceAlpha1=0.25; 
edgeWidth=2; 
markerSize=50; 
cMap=gjet(6);

%% Examples 
% 

%% Converting a hexahedral element to tetrahedral elements

%%
% Creating an example hexahedral element
V=[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1;]; %nodes
E=1:8; %Element

%%

[F,CF]=element2patch(E);  %Patch data for plotting

cFigure;
subplot(2,3,1); hold on;
title('Original element set','FontSize',fontSize);
gpatch(F,V,cMap(1,:),'k',1);
patchNormPlot(F,V,0.25);
colormap(cMap);
axisGeom(gca,fontSize);

for q=1:1:5
    
    % Subdeviding the hexahedral element
    convertMethod=q; %Corresponse 
    [Es,Vs,Cs]=hex2tet(E,V,[],convertMethod);
    [Fs,CFs]=element2patch(Es); %Patch data for plotting
    
    subplot(2,3,q+1); hold on;
    title(['Converted, method: ',num2str(q)],'FontSize',fontSize);
    gpatch(Fs,Vs,cMap(q+1,:),'k',0.5);
    patchNormPlot(Fs,Vs,0.25);
    
    colormap(cMap);
    axisGeom(gca,fontSize);
end
drawnow;

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
