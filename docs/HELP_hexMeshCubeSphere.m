clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=15;
faceAlpha1=0.5;

%%

cPar.boxWidth=1;
cPar.OuterSphereRadius=cPar.boxWidth/6; 
cPar.InnerSphereRadius=cPar.OuterSphereRadius/2; 
cPar.CoreSphereRadius=cPar.InnerSphereRadius/2; 
cPar.numElementsCube=16;
cPar.numElementsCubeSphere=6;
cPar.numElementsSphereMantel=6;
cPar.numElementsSphereCore=6;
cPar.nSmooth=15;

[E,V,CE,Fb,Cb]=hexMeshCubeSphere(cPar);

%%

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),CE(L,:),'hex8');

%%

cFigure;
subplot(1,2,1); hold on;
title('Cut-view of the mesh','FontSize',fontSize);

gpatch(Fs,V,Cs);

axisGeom(gca,fontSize);
colormap(gca,gjet(3)); icolorbar; 
camlight headlight;

subplot(1,2,2); hold on;
title('Mesh boundaries','FontSize',fontSize);

gpatch(Fb,V,Cb,'none',0.5);

axisGeom(gca,fontSize);
colormap(gca,gjet(8)); icolorbar; 
camlight headlight;

drawnow; 

