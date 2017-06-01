%% tet4_tet10
% Below is a demonstration of the features of the |tet4_tet10| function

%%
clear; close all; clc;

% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=0.5;
edgeColor='k';
edgeWidth1=2;
edgeWidth2=1;
markerSize1=75;
markerSize2=10;

%% CONVERSION FROM TET4 TO TET10, EXAMPLE FOR A SINGLE TETRAHEDRON
% Creating a single 4-node tetrahedron
[V4,~]=platonic_solid(1,1); %q indicates solid type, r is the radius
TET4=[1:4];
[F4,~]=element2patch(TET4,[],'tet4');

%%
% Converting to a single 10-node tetrahedron
[TET10,V10,~]=tet4_tet10(TET4,V4,{});
[F10,~]=element2patch(TET10,[],'tet10');

%%
% Plotting elements
hf=cFigure; % Open figure for plotting
subplot(1,2,1); hold on;
title('A linear tetrahedron','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F4,'Vertices',V4);
set(hp,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth1,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize1);
[hn]=patchNormPlot(F4,V4,0.25); %Plotting face normals

for q=1:1:size(TET4,2)
    text(V4(TET4(1,q),1),V4(TET4(1,q),2),V4(TET4(1,q),3),[' ',num2str(q)],'FontSize',50);
end

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;

subplot(1,2,2); hold on;
title('A quadratic tetrahedron','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F10,'Vertices',V10);
set(hp,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth1,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize1);
[hn]=patchNormPlot(F10,V10,0.25); %Plotting face normals

for q=1:1:size(TET10,2)
    text(V10(TET10(1,q),1),V10(TET10(1,q),2),V10(TET10(1,q),3),[' ',num2str(q)],'FontSize',50);
end

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;

drawnow; 

%% CONVERSION FROM TET4 TO TET10, EXAMPLE FOR A TETRAHEDRON MESH WITH NODAL PARAMETERS
n=3;
for q=1:1:n
    [TET4,V4]=subTet(TET4,V4,1);
end

[F4,~]=element2patch(TET4,[],'tet4');

%Create nodal result e.g. displacement and color
V4d=V4;
V4d(:,1)=V4d(:,1)+0.5.*sin(pi*V4d(:,3));
V4d(:,2)=V4d(:,2)+0.5.*cos(pi*V4d(:,3));
D4=V4d-V4;

C4=V4(:,1); %Color towards X

%%
% Converting to a single 10-node tetrahedron
tet4_cell={D4,C4};
[TET10,V10,tet10_cell]=tet4_tet10(TET4,V4,tet4_cell);
[F10,~]=element2patch(TET10,[],'tet10');
D10=tet10_cell{1};
C10=tet10_cell{2};
V10d=V10+D10;

%%
% Plotting elements
hf=cFigure; % Open figure for plotting
subplot(2,2,1); hold on;
title('A linear tetrahedron mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F4,'Vertices',V4);
set(hp,'FaceColor','g','FaceAlpha',faceAlpha,'lineWidth',edgeWidth2,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize2);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,3); hold on;
title('Data on linear tetrahedron mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F4,'Vertices',V4d);
set(hp,'FaceColor','flat','CData',C4,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth2,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize2);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,2); hold on;
title('A quadratic tetrahedron mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F10,'Vertices',V10);
set(hp,'FaceColor','r','FaceAlpha',faceAlpha,'lineWidth',edgeWidth2,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize2);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,4); hold on;
title('Mapped data on quadratic tetrahedron mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

hp=patch('Faces',F10,'Vertices',V10d);
set(hp,'FaceColor','flat','CData',C10,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth2,'edgeColor',edgeColor,'Marker','.','MarkerSize',markerSize2);
colormap jet; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
view([-50,12])
camlight('headlight'); lighting flat;

drawnow; 

%%


%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
