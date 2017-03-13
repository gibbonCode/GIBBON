%% tet2hex
% Below is a demonstration of the features of the |tet2hex| function

%% Syntax
% |[Es,Vs]=tet2hex(E,V);|

%% Description
%
%% Examples

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceColor1='g';
faceColor2='r';
faceAlpha1=0.3;
faceAlpha2=1;
edgeColor=0.*ones(1,3);
edgeWidth=2;
markerSize=2;
cMap=gjet(250);

%% Example converting a single tetrahedron to 4 hexahedrons

%%
% Creating an example tetrahedron
[V,~]=platonic_solid(1,1);
E=[1:4];

%%
% Convert tetrahedron to hexahedral elements

[Es,Vs]=tet2hex(E,V);

%% Visualization

[F]=element2patch(E);  %Patch data for plotting

Cs=(1:1:size(Es,1))';
[Fs,CFs]=element2patch(Es,Cs); %Patch data for plotting

cFigure;
subplot(1,2,1); 
title('Original tetrahedral element','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);
view(3); grid on; axis equal; axis tight;
set(gca,'FontSize',fontSize);

subplot(1,2,2); 
title('Converted hexahedral elements','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',Fs,'Vertices',Vs,'FaceColor','flat','CData',CFs,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);

colormap(cMap);
view(3); grid on; axis equal; axis tight;
set(gca,'FontSize',fontSize);

drawnow;

%% Example converting a set of tetrahedral elements

%%
% Creating an example set of hexahedrons

[V,~]=platonic_solid(1,1);
E=[1 2 4 3];

n=0; 
if n>0
    for q=1:1:n
        [E,V]=subTet(E,V,1);
    end
end

C=(1:1:size(E,1))';

%%
% Subdeviding the hexahedral element
[Es,Vs]=tet2hex(E,V);
 
%% Visualization

[F,CF]=element2patch(E,C);  %Patch data for plotting

[Fs,CFs]=element2patch(Es,repmat(C,[4 1])); %Patch data for plotting

cFigure;
subplot(1,2,1); 
title('Original tetrahedral element set','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);
hp=patchNormPlot(F,V);
view(3); grid on; axis equal; axis tight;
set(gca,'FontSize',fontSize);

subplot(1,2,2); 
title('Converted hexahedral elements','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',Fs,'Vertices',Vs,'FaceColor','flat','CData',CFs,'EdgeColor',edgeColor,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth);
hp=patchNormPlot(Fs,Vs);
colormap(cMap);
view(3); grid on; axis equal; axis tight;
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