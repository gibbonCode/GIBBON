%% voronoiDiagramEven2D
% Below is a basic demonstration of the features of the |voronoiDiagramEven2D| function.

%% Syntax
% |[Vv,Fv]=voronoiDiagramEven2D(DT,numPointsVoronoi);|

%% Description
% Samples the cells of a Voronoi tesselation using the same amount of
% points

%% Examples

%%
clear; close all; clc;

% PLOT SETTINGS
markerSize=15;
lineWidth=2;
fontSize=15; 
fAlpha1=0.25; 
fAlpha2=0.8; 
lineWidth1=1; 
lineWidth2=3; 

%% Example: Sampling Voronoi cells with the same amount of points

%%
% EXAMPLE DELAUNAY TRIANGULATION

%Boundary and mesh parameters
ns=50; %Number of points on outer boundary (defines how well the circle is sampled)
rOut=1; %Outer radius of circular boundary
pointSpacing=rOut/5; %Approximate initial point spacing for point seeding
stdP=pointSpacing/2*ones(1,2); %Standard deviations for random point offset after point seeding

%Creating boundary curve
tt=linspace(0,2*pi,ns);
tt=tt(1:end-1);
r=rOut.*ones(size(tt));
[x,y] = pol2cart(tt,r);
Vb=[x(:) y(:)];

%Create Delaunay derived mesh
regionCell={Vb};
[Ft,Vt,~,DT]=regionTriMeshRand2D(regionCell,pointSpacing,stdP,1,0);

%% 
% COMPUTE NORMAL VORONOI TESSELATION REQUIRING LOOP FOR PLOTTING

[Vv,Fv_cell] =voronoiDiagram(DT);
Vv(isinf(Vv))=NaN;

hf1=cFigure;
title('Normal Voronoi tesselation','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
hold on;
patch('faces',Ft,'vertices',Vt,'FaceColor','r','FaceAlpha',fAlpha1,'LineWidth',lineWidth1);

HP=cellPatch(Vv,Fv_cell,'g');
set(HP,'LineWidth',lineWidth,'faceAlpha',fAlpha2,'LineWidth',lineWidth2);

% for q=1:1:numel(Fv_cell)
%     fv=Fv_cell{q};
%     vv=Vv(fv,:);
%     if all(~isnan(vv))
%         patch('faces',fv,'vertices',Vv,'FaceColor','g','faceAlpha',fAlpha2,'LineWidth',lineWidth2);      
%     end
% end

axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;
drawnow;

%%
% Use |voronoiDiagramEven2D| to sample all cells with same amount of points

numPointsVoronoi1=3; 
numPointsVoronoi2=5; 
numPointsVoronoi3=11; 
numPointsVoronoi4=25; 

[Vv1,Fv1]=voronoiDiagramEven2D(DT,numPointsVoronoi1);
[Vv2,Fv2]=voronoiDiagramEven2D(DT,numPointsVoronoi2);
[Vv3,Fv3]=voronoiDiagramEven2D(DT,numPointsVoronoi3);
[Vv4,Fv4]=voronoiDiagramEven2D(DT,numPointsVoronoi4);

%%
% Now for loops can be avoided for plotting

hf1=cFigure;
subplot(2,2,1);
title(['Sampled using ',num2str(numPointsVoronoi1),' points'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);
hold on;
patch('faces',Ft,'vertices',Vt,'FaceColor','r','FaceAlpha',fAlpha1,'LineWidth',lineWidth1);
patch('faces',Fv1,'vertices',Vv1,'FaceColor','g','FaceAlpha',fAlpha2,'LineWidth',lineWidth2);
axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;

subplot(2,2,2);
title(['Sampled using ',num2str(numPointsVoronoi2),' points'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);
hold on;
patch('faces',Ft,'vertices',Vt,'FaceColor','r','FaceAlpha',fAlpha1,'LineWidth',lineWidth1);
patch('faces',Fv2,'vertices',Vv2,'FaceColor','g','FaceAlpha',fAlpha2,'LineWidth',lineWidth2);
axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;

subplot(2,2,3);
title(['Sampled using ',num2str(numPointsVoronoi3),' points'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);
hold on;
patch('faces',Ft,'vertices',Vt,'FaceColor','r','FaceAlpha',fAlpha1,'LineWidth',lineWidth1);
patch('faces',Fv3,'vertices',Vv3,'FaceColor','g','FaceAlpha',fAlpha2,'LineWidth',lineWidth2);
axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;

subplot(2,2,4);
title(['Sampled using ',num2str(numPointsVoronoi4),' points'],'FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);
hold on;
patch('faces',Ft,'vertices',Vt,'FaceColor','r','FaceAlpha',fAlpha1,'LineWidth',lineWidth1);
patch('faces',Fv4,'vertices',Vv4,'FaceColor','g','FaceAlpha',fAlpha2,'LineWidth',lineWidth2);
axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
