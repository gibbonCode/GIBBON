%% foamWrap
% Below is a demonstration of the features of the |foamWrap| function

%% Syntax
% |[FT,VT,CT,CT_c]=foamWrap(F,V,C,cPar);|

%% Description
% Use |foamWrap| to generate a foam like structure on top of an input mesh

%% Examples

clear; close all; clc;

%% 
% Plot Settings

fontSize=15;
faceAlpha=1;
edgeColor=0.1*ones(1,3);
edgeWidth=1;
cmap=gjet(250);

%% 
% Create surface model 

[F,V,~]=geoSphere(0,1); %Geodesic sphere
% [F,V]=parasaurolophus;
% [F,V]=cow;
% [F,V]=graphicsModels(4);
% [F,V]=stanford_bunny;

[F,V,C,indIni]=triPolyDualRefine(F,V);

%%
cFigure;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','k');

camlight headlight; 
colormap(cmap);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on; axis off; 


%%

cPar.n=3; 
cPar.dirFlip=1; 
cPar.foamThickness=[]; %Empty uses default which is mean edgelength based
cParSmooth.Method='HC';
cParSmooth.n=25;
cPar.cParSmooth=cParSmooth; 

%%
L_remove=true(size(F,1),1);
[FT,VT,CT,CT_c]=foamWrap(F,V,C,cPar);

%%
cFigure; hold on; 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',FT(CT_c==1,:),'Vertices',VT,'FaceColor','flat','CData',CT_c(CT_c==1),'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor','none');
hp=patch('Faces',FT(CT_c~=1,:),'Vertices',VT,'FaceColor','flat','CData',CT_c(CT_c~=1),'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor','none');

camlight headlight; 
% lighting phong;
colormap(gjet(numel(unique(CT_c)))); colorbar;


set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
axis off;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>