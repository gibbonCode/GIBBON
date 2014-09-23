%% minDist
% Below is a demonstration of the features of the |minDist| function

%%
closeFall;

%%
% PLOT SETTINGS
fig_color='w'; fig_colordef='white';
font_size=15;
cmap=gray(250);
falpha=1;
patch_types={'sx','sy','sz','v'};
ptype=3;
no_slices=4;
mark_siz1=20;
mark_siz2=5;
mark_siz3=50;
line_width1=2;
F_alpha1=0.5;
F_alpha2=1;

%% EXAMPLE FOR POINT CLOUD OR SURFACE DISTANCE COMPUTATION

%% 
% Building test surfaces

%Defining shape 1 as a sphere
[F1,V1,~]=geoSphere(2,1);

%Defining shape 2 as a deformed sphere
[F2,V2,Vs]=geoSphere(3,1);
freqDef=3;
ampDef=0.25;
ampDefDiff=0.25;
n1=Vs(:,3)+(ampDef-ampDefDiff)+ampDef*sin(freqDef*Vs(:,1));
[V2(:,1),V2(:,2),~]=sph2cart(Vs(:,1),Vs(:,2),n1);

%%
% Plotting surfaces

hf1=figuremax(fig_color,fig_colordef);
title('The two surfaces','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F1,'vertices',V1,'FaceColor','g','FaceAlpha',F_alpha1);
patch('faces',F2,'vertices',V2,'FaceColor','b','FaceAlpha',F_alpha1);

axis equal; view(3); axis tight; grid on;
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;

%% 
% Get closest point based distance metric
D2=minDist(V2,V1);

%%
% On this type of use see also the |triSurfSetDist| function

%%
% Plotting results

[CF]=vertexToFaceMeasure(F2,D2);

hf2=figuremax(fig_color,fig_colordef);
title('Closest point distance metric on surface 2','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2,'vertices',V2,'FaceColor','flat','CData',CF);
patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',F_alpha1,'EdgeColor','None');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;

%% EXAMPLE FOR NEAREST NEIGHBOUR SEARCH AND INTERPOLATION

%%
% Building test surfaces

%Defining shape 1 as a sphere
[F1,V1,~]=geoSphere(3,1);

%Defining shape 2 as a denser sphere
[F2,V2,Vs]=geoSphere(5,1);

%Simulate some kind of result on the coarse sphere
C1=triplyPeriodicMinimal(6.*V1,'g');

%%
% Find nearest neighbours
[~,indMin]=minDist(V2,V1);

%%
% Interpolation now reduces to simple indexing into the array
C2=C1(indMin);

%%
% Get "true" color to compare
C2_true=triplyPeriodicMinimal(6.*V2,'g');

%%
% Plotting surfaces
[CF1]=vertexToFaceMeasure(F1,C1);
[CF2]=vertexToFaceMeasure(F2,C2);
[CF2_true]=vertexToFaceMeasure(F2,C2_true);

hf2=figuremax(fig_color,fig_colordef);
subplot(1,3,1);
title('Coarse input','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F1,'vertices',V1,'FaceColor','flat','CData',CF1,'FaceAlpha',F_alpha2,'EdgeColor','k');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 

subplot(1,3,2);
title('Upsampled using nearest neighbours','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2,'vertices',V2,'FaceColor','flat','CData',CF2,'FaceAlpha',F_alpha2,'EdgeColor','k');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 

subplot(1,3,3);
title('Truth','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2,'vertices',V2,'FaceColor','flat','CData',CF2_true,'FaceAlpha',F_alpha2,'EdgeColor','none');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 

drawnow;

%%
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)
