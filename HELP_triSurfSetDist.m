%% triSurfSetDist
% Below is a demonstration of the features of the |triSurfSetDist| function

%%
closeFall;

%%
% PLOT SETTINGS
fig_color='w'; fig_colordef='white';
font_size=20;
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

%% CLOSEST POINT BASED DISTANCE METRIC

[D2]=triSurfSetDist(F2,V2,F1,V1,'dist');

%%
% The above is equivalent to: 
D2=minDist(V2,V1);

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

%% RAY TRACING DISTANCE METRIC

[D2]=triSurfSetDist(F2,V2,F1,V1,'ray');

%%
% Plotting results

[CF]=vertexToFaceMeasure(F2,D2);
L=~isnan(CF); %Check for NaN's

hf3=figuremax(fig_color,fig_colordef);
title('Ray-traced distance metric on surface 2','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2(L,:),'vertices',V2,'FaceColor','flat','CData',CF(L));
patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',F_alpha1,'EdgeColor','None');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;

%% NOTE ON RAY TRACING DISTANCE METRIC AND NORMAL DIRECTIONS THAT DO NOT INTERSECT
% The 'dist' method always finds a closest distance result for all points
% of the surface. However for the ray method the surface normals do not
% always trace to the other surface. The below example illustrates this.
% The curvature of the top and bottom regions means the normal direction
% rays do not intersect with the surface 1. 

%Stretching shape in Z-direction 
V2(:,3)=V2(:,3)*2;

%Compute distance metric
[D2]=triSurfSetDist(F2,V2,F1,V1,'ray');

%%
% Plotting results

[CF]=vertexToFaceMeasure(F2,D2);
L=~isnan(CF); %Check for NaN's

hf3=figuremax(fig_color,fig_colordef);
title('Ray-traced distance metric on surface 2','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2(L,:),'vertices',V2,'FaceColor','flat','CData',CF(L));
patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',F_alpha1,'EdgeColor','None');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;

%% EXAMPLE TO FIX NAN's

%Logic for NaN's
L=isnan(D2); 

%Use dist method where NaN occured
[D2_nan]=triSurfSetDist(F2,V2(L,:),F1,V1,'dist');
D2(L)=D2_nan; 

%%
% Plotting results

[CF]=vertexToFaceMeasure(F2,D2);

hf3=figuremax(fig_color,fig_colordef);
title('Combined distance metric on surface 2','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2,'vertices',V2,'FaceColor','flat','CData',CF);
patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',F_alpha1,'EdgeColor','None');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;

%% EXAMPLE OF "HYBRID APPROACH" 
% A hybrid approach is also possible whereby the output is the smallest
% distance for the two methods (nanmin is used so NaN values due to ray
% tracing enforces the distance method on these locations instead). 

%Compute distance metric
[D2]=triSurfSetDist(F2,V2,F1,V1,'dist-ray');

%%
% Plotting results
[CF]=vertexToFaceMeasure(F2,D2);

hf3=figuremax(fig_color,fig_colordef);
title('Minimum metric between "closest point" and "ray-traced distance" on surface 2','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size);zlabel('Z','FontSize',font_size); 
hold on; 
patch('faces',F2,'vertices',V2,'FaceColor','flat','CData',CF);
patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',F_alpha1,'EdgeColor','None');

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off; 
set(gca,'FontSize',font_size); 
camlight headlight; 
drawnow;




