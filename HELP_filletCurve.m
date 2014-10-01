%% filletCurve
% Below is a demonstration of the features of the |filletCurve| function
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

%% FILLETING A CURVE IN 3D

%%
% Simulating a curve with sharp features
Vt=[0 0 0; 10 0 0; 5 10 0; 10 0 10; 0 10 10; ];

%%
%Setting control parameters
r=2; %Fillet radius
np=25; %Number of points used to construct each fillet edge
closedLoopOption=0; %Use 1 if curve represents a closed loop but containes unique points
[VN]=filletCurve(Vt,r,np,closedLoopOption);

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);
title('A filleted curve','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Vt,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(VN,'r.-','lineWidth',lineWidth2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% FILLETING A CLOSED CURVE IN 3D
closedLoopOption=1; 
[VN]=filletCurve(Vt,r,np,closedLoopOption);

%%
% Plotting results
hf2=figuremax(figColor,figColorDef);
title('A filleted curve based on closed end conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Vt,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV([Vt(1,:);Vt(end,:)],'g.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(VN,'r.-','lineWidth',lineWidth2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% EXTRUDING A FILLETED CURVE FOR CAD LIKE MODEL BUILDING

%Sketching side profile
x=[0 0 5 15 15];
y=[0 4 9 10 0];
V=10*[x(:) y(:)];

%Fillet sketch
r=15; %Fillet radius
np=50; %Number of points used to construct each fillet edge
closedLoopOption=0; %Use 1 if curve represents a closed loop but containes unique points
[Vc]=filletCurve(V,r,np,closedLoopOption);

%%
% Plotting sketch
hf1=figuremax(figColor,figColorDef);
title('The side profile sketch','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(V,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(Vc,'r-','lineWidth',lineWidth2,'MarkerSize',markerSize1);

axis equal; view(2); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%%
% Extruding model
cPar.pointSpacing=10;
cPar.depth=450; 
cPar.patchType='tri'; 
cPar.dir=0;
[F_tri,V_tri]=polyExtrude(Vc,cPar);

%% 
% Plotting meshed model
hf2=figuremax(figColor,figColorDef);
title('The extruded model mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

hp=patch('faces',F_tri,'Vertices',V_tri);
% [hp2]=patchNormPlot(F_tri,V_tri,2*pointSpacing);

set(hp,'FaceColor','g','EdgeColor','k','FaceAlpha',faceAlpha,'LineWidth',lineWidth3);
camlight headlight;
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>