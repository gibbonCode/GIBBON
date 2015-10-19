%% patch_dual
% Below is a demonstration of the features of the |patch_dual| function

%%
clear; close all; clc; 

%%
% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=0.5;
plotColor1=0.2.*ones(1,3);
plotColor2=0.5.*ones(1,3);
edgeWidth=3;
markerSize=15;
cmap=autumn(250);
figStruct.ColorDef='black';
figStruct.Color='k';

%% EXAMPLE: The "Buckminster Fuller" dome triangulation and its dual
% The patch_dual function assumes that a valid and appropriate dual exists
% for the input patch data specified by F and V (faces and vertices). If
% they are not appropriate the output may for instance not form an
% enclosing shape or output faces may not be planar. 

%Defining geodesic dome triangulation
r=1; %sphere radius
n=2; %Refinements   
[F,V,~]=geoSphere(n,r);

%Deriving the dual of the patch
[Vd,Fd]=patch_dual(V,F);

%Plotting results
hf=cFigure(figStruct);
hold on;
% title('A geodesic sphere triangulation and its dual consisting of pentagons and hexagons','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','none','EdgeColor',plotColor2,'LineWidth',edgeWidth,'Marker','o','MarkerFaceColor',plotColor1,'MarkerEdgeColor','none','MarkerSize',markerSize);

%Splitting up the plotting due to difference is face types (e.g.
%pentagons,or hexagons)
for i=1:1:numel(Fd);
    Ft=Fd{i};
    hp=patch('Faces',Ft,'Vertices',Vd);
    C=rand(size(Ft,1),1); %Create random color
    set(hp,'FaceColor','flat','CData',C,'FaceAlpha',0.6,'EdgeColor',plotColor2,'LineWidth',edgeWidth,'Marker','o','MarkerFaceColor',plotColor2,'MarkerEdgeColor','none','MarkerSize',markerSize);
end
colormap(cmap);
axis equal; axis tight; view(3); axis vis3d; axis off; 
set(gca,'FontSize',fontSize);
camlight headlight; lighting flat; 
drawnow;

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
