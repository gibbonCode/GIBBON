%% geoSphere
% Below is a demonstration of the features of the |geoSphere| function

%%
clear all; close all; clc; 

%% 
% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceAlpha=0.75;
edgeColor=0.3*ones(1,3);
edgeWidth=1.5;

%% Building a geodesic dome
% The function inputs are n and r which define the mesh refinement and
% radius respectively. The mesh refinement number n defines the number of
% subtriangulation (see function |subTri|) iterations performed on an
% icosahedron. Below is a visualisation for n=0:1:3. The function outputs
% the geodesic dome faces (F) and vertices (V) and also the spherical
% coordinates of the vertices (Vs) (this output is suppressed in the
% example below).

% Open figure for plotting
hf=figuremax(fig_color,fig_colordef); 

%Defining geodesic dome
r=1; %sphere radius
n=0:1:3; %Refinements   
pColors=autumn(numel(n));
for q=1:1:numel(n);
    [F,V,~]=geoSphere(n(q),r); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' refinement iterations'],'FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hp=patch('Faces',F,'Vertices',V);
    set(hp,'FaceColor',pColors(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    camlight headlight; 
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  grid on;
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)