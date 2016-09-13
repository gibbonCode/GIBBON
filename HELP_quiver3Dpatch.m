%% quiver3Dpatch
% Below is a demonstration of the features of the |quiver3Dpatch| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white'; 
cMap=jet(250);
faceAlpha1=1;
faceAlpha2=1;
edgeColor1='none';
edgeColor2='none';
cMap1=jet(250);
cMap2=gray(250);
fontSize=8; 

%% Plotting a vector
% Below is a visualisation of the basec vector style

%Defining a single vector colinear with the Z-axis with length 2
X=0; Y=0; Z=0; %Vector origin (position vector components)
u=0; v=0; w=2; %Vector components
G=sqrt(u.^2+v.^2+w.^2); %Vector magnitude
cLim=[0 max(G(:))];
Cv=[]; %If empty then vector magnitude based scaling is used

a=[min(G(:)) max(G(:))]; %Arrow length scaling to magnitude range
[F1,V1,C1]=quiver3Dpatch(X,Y,Z,u,v,w,Cv,a);

h1=figuremax(fig_color,fig_colordef);
title('Basic vector style using 7 vertices and 6 faces');
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
patch('Faces',F1,'Vertices',V1,'EdgeColor','k', 'CData',C1,'FaceColor','flat','FaceAlpha',0.5,'Marker','.','MarkerSize',25); 
colormap(cMap1); colorbar; caxis(cLim);
view(3); grid on; axis equal; axis tight; axis vis3d; 
set(gca,'FontSize',fontSize);

%% Defining vector lengths and colours

a=[min(G(:)) max(G(:))]; %Arrow length scaling to magnitude range
[F1,V1,C1]=quiver3Dpatch(X,Y,Z,u,v,w,Cv,a);

a=[1 1]; %Arrow length scaling min=1, max=1
[F2,V2,C2]=quiver3Dpatch(X,Y,Z,u,v,w,Cv,a);

a=[min(G(:)) max(G(:))]; %Arrow length scaling to magnitude range
Cv=zeros(size(X));
[F3,V3,C3]=quiver3Dpatch(X,Y,Z,u,v,w,Cv,a);

C4=gray2RGBColorMap(C3,cMap2,cLim);

h1=figuremax(fig_color,fig_colordef);
subplot(2,2,1);
title('Vector with length and color according to magnitude');
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
patch('Faces',F1,'Vertices',V1,'EdgeColor','k', 'CData',C1,'FaceColor','flat','FaceAlpha',1); 
colormap(cMap1); colorbar; caxis(cLim);
view(3); grid on; axis equal; axis tight; axis vis3d; 
set(gca,'FontSize',fontSize);

subplot(2,2,2);
title('Vector with a scaled length but color according to magnitude');
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
patch('Faces',F2,'Vertices',V2,'EdgeColor','k', 'CData',C2,'FaceColor','flat','FaceAlpha',1); 
colormap(cMap1); colorbar; caxis(cLim);
view(3); grid on; axis equal; axis tight; axis vis3d; 
set(gca,'FontSize',fontSize);

subplot(2,2,3);
title('Vector with length according to magnitude a user specified colormapo driven color');
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
patch('Faces',F3,'Vertices',V3,'EdgeColor','k', 'CData',C3,'FaceColor','flat','FaceAlpha',1); 
colormap(cMap1); colorbar; caxis(cLim);
view(3); grid on; axis equal; axis tight; axis vis3d; 
set(gca,'FontSize',fontSize);

subplot(2,2,4);
title('Vector with length according to magnitude and user specified RGB driven color');
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
patch('Faces',F3,'Vertices',V3,'EdgeColor','k', 'FaceVertexCData',C4,'FaceColor','flat','FaceAlpha',1); 
view(3); grid on; axis equal; axis tight; axis vis3d; 
set(gca,'FontSize',fontSize);
camlight headlight; lighting phong
drawnow;

%% Example visualising coordinate system base vectors

originBasis1=[0 0 0];
E1=eye(3,3);
C1=[2 1 0];

originBasis2=[0 0 0];
E2=[2/3 -1/3  2/3; 2/3 2/3 -1/3; -1/3 2/3 2/3];
C2=[5 4 3];

[Fc1,Vc1,Cc1]=quiver3Dpatch(originBasis1(1)*ones(1,3), originBasis1(2)*ones(1,3), originBasis1(3)*ones(1,3),E1(:,1),E1(:,2),E1(:,3),C1',[1 1]);
[Fc2,Vc2,Cc2]=quiver3Dpatch(originBasis2(1)*ones(1,3), originBasis2(2)*ones(1,3), originBasis2(3)*ones(1,3),E2(:,1),E2(:,2),E2(:,3),C2',[1 1]);

h1=figuremax(fig_color,fig_colordef);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
title('Visualizing base vectors','FontSize',fontSize);
hp1=patch('Faces',Fc1,'Vertices',Vc1,'EdgeColor','k','FaceColor','flat','FaceVertexCData',Cc1,'FaceAlpha',1); hold on;
hp1=patch('Faces',Fc2,'Vertices',Vc2,'EdgeColor','k','FaceColor','flat','FaceVertexCData',Cc2,'FaceAlpha',0.5); hold on;
view(3); grid on; axis equal; axis vis3d; view([137.5,24]);
set(gca,'FontSize',fontSize);
colormap jet; 
drawnow;

%% Example visualising face normals of patch data
hf=figuremax(fig_color,fig_colordef); % Open figure for plotting

[Fs,Vs,~]=geoSphere(2,1);
title('Displaying face normals','FontSize',fontSize);
hp=patch('Faces',Fs,'Vertices',Vs,'FaceColor','g');

%Plotting face normals
[hn]=patchNormPlot(Fs,Vs,0.3);

set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting flat;

%% Example for multidimensional image data 1: colormap driven vectors combined with RGB driven iso-surfaces

%%
% Simulating 3D volume and vector data
n=27;
[X,Y,Z]=meshgrid(linspace(-4.77,4.77,n));
phi=(1+sqrt(5))/2;
M=2 - (cos(X + phi*Y) + cos(X - phi*Y) + cos(Y + phi*Z) + cos(Y - phi*Z) + cos(Z - phi*X) + cos(Z + phi*X));

%%
% Simulating vector data 
%Vector data here based on the gradient of the image
[u,v,w] = gradient(M); 
G=hypot(hypot(u,v),w); %Vector lenghts

%Iso-surface patch data to illustrate joint plotting
c_iso1=0; c_iso2=5;
[Fi1,Vi1,Ci1] = isosurface(X,Y,Z,M,c_iso1,M); 
[Fi2,Vi2,Ci2] = isosurface(X,Y,Z,M,c_iso2,M); 
a=[min(G(:)) max(G(:))]; %Arrow length scaling
L=G>0.9; %Logic indices for arrows
[Fv,Vv,Cv]=quiver3Dpatch(X(L),Y(L),Z(L),u(L),v(L),w(L),G(L),a);

cLim=[min(M(:)) max(M(:))]; %Colorbar limits
[Ci1n]=gray2RGBColorMap(Ci1,cMap2,cLim);
[Ci2n]=gray2RGBColorMap(Ci2,cMap2,cLim);

h1=figuremax(fig_color,fig_colordef);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
title('Colormap driven vector colors and RGB driven isosurfaces','FontSize',fontSize);
patch('Faces',Fv,'Vertices',Vv,'EdgeColor',edgeColor1, 'CData',Cv,'FaceColor','flat','FaceAlpha',1); 
patch('Faces',Fi1,'Vertices',Vi1,'FaceColor','flat','FaceVertexCData',Ci1n,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2); hold on;
patch('Faces',Fi2,'Vertices',Vi2,'FaceColor','flat','FaceVertexCData',Ci2n,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2); hold on;
colormap(cMap1); colorbar; caxis([min(Cv(:)) max(Cv(:))]);
view(3); grid on; axis equal; axis vis3d; 
set(gca,'FontSize',fontSize);
camlight headlight; lighting phong
drawnow;

%% Example for multidimensional image data 2: RGB driven vectors combined with colormap driven iso-surfaces
% Angle driven color can also be specified e.g. RGB values indicating vector angle

%Specifying angle dependant RGB type color
Xc=repmat(u(L),[6,1]); Yc=repmat(v(L),[6,1]); Zc=repmat(w(L),[6,1]);
Crgb=[Xc(:) Yc(:) Zc(:)];
M=sqrt(Crgb(:,1).^2+Crgb(:,2).^2+Crgb(:,3).^2);
Crgb=abs(Crgb./(M*ones(1,3))); %Normalising color

%%
% Defining a sphere to show the color mapping
[F,V,~]=geoSphere(4,1);
Xs=V(:,1); Ys=V(:,2); Zs=V(:,3);
C=[mean(Xs(F),2) mean(Ys(F),2) mean(Zs(F),2)]; %color for angles

%%
% The figure now demonstrates isosurfaces for the image data with overlain
% the gradient vectors coloured according to their direction

h1=figuremax(fig_color,fig_colordef);
subplot(1,2,1);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
title('RGB driven vector colors and colormap driven isosurfaces','FontSize',fontSize);

Cv=vecnormalize(Vv);
patch('Faces',Fv,'Vertices',Vv,'EdgeColor','none', 'FaceVertexCData',Crgb,'FaceColor','flat','FaceAlpha',1); 
patch('Faces',Fi1,'Vertices',Vi1,'FaceColor','flat','CData',Ci1,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2); hold on;
patch('Faces',Fi2,'Vertices',Vi2,'FaceColor','flat','CData',Ci2,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2); hold on;
view(3); grid on; axis equal; axis vis3d; 
set(gca,'FontSize',fontSize);
colormap(cMap2); colorbar; 
camlight headlight; lighting phong

subplot(1,2,2);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','flat','FaceVertexCData',abs(C),'EdgeColor','none','FaceAlpha',1);
DCM=eye(3,3);
origin=[0 0 0];
[Fa,Va,Ca]=quiver3Dpatch(origin(1)*ones(1,3), origin(2)*ones(1,3), origin(3)*ones(1,3),-DCM(:,1),-DCM(:,2),DCM(:,3),[],[3,3]);
hp2=patch('Faces',Fa,'Vertices',Va,'EdgeColor','k','FaceColor','flat','FaceVertexCData',repmat(eye(3,3),6,1),'FaceAlpha',1); hold on;

view(3); axis tight; axis square; axis vis3d; view(-45,30);
set(gca,'FontSize',fontSize); drawnow;
axis off; 
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
