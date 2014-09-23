%% tet4_tet10
% Below is a demonstration of the features of the |tet4_tet10| function

%%
clear all; close all; clc;

% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=0.75;
edgeColor='k';
edgeWidth=2;
markerSize=5;

%%




%Defining the faces (F) and vertices (V) of a platonic solid
[V,~]=platonic_solid(1,1); %q indicates solid type, r is the radius
TET4=[1 2 3 4];
[F,~]=element2patch(TET4,[]);

hf=figuremax(fig_color,fig_colordef); % Open figure for plotting
subplot(1,2,1); hold on;
title('A linear tetrahedron','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);

%Plotting face normals
[hn]=patchNormPlot(F,V,0.5);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;


%%
%
% <<gibbVerySmall.gif>>
%
% GIBBON
%
% Kevin M. Moerman (kevinmoerman@hotmail.com)