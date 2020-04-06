clear; close all; clc;

%% Create example geometry

r=1;
[F,V]=hemiSphereMesh(2,r,0);
Eb=patchBoundary(F,V);
VF=patchCentre(F,V);
logicKeep=VF(:,3)<r/2;
F=F(logicKeep,:);
[F,V,indFix]=patchCleanUnused(F,V);
Eb=indFix(Eb);

numEdgesNeeded=size(Eb,1)+5; 

%% 
% Split edges until number of desired edges is reached

[Fn,Vn,Ebn]=triSurfSplitBoundary(F,V,Eb,numEdgesNeeded);

%%
% Visualize result

cFigure; 
subplot(1,2,1);
gpatch(F,V,'bw','k');
gpatch(Eb,V,'none','g',1,3);
axisGeom;
camlight headlight; 

subplot(1,2,2);
gpatch(Fn,Vn,'gw','k');
gpatch(Ebn,Vn,'none','b',1,3);
axisGeom;
camlight headlight; 

drawnow; 


