%% linspacen
% Below is a demonstration of the features of the |linspacen| function

%% Summary 
% This function is a generalization of the |linspace| function to N
% dimensions. The output C is a matrix of size [size(A) n] such that "it
% goes" from A to B in n steps in the last dimention. The input variables A
% and B (scalars, vectors or matrices). For scalar input this function is
% equivalent to |linspace| (but slower due to repmat operation). 
% Clearly the inputs A and B should have the same size.
%
% See also: |linspace| and |subCurve|
%%

clear; close all; clc;

%%
% PLOT SETTINGS
figColor='w'; 
figColorDef='white';
fontSize=10;
cMap=jet(250);

%% Control parameters for examples
% Number of steps used in examples
n=6; 
%% Use of |linspacen| for scalars
% For scalar input |linspacen| is equivalent to |linspace|

A=0 
B=1 
C=linspacen(A,B,n) 

%% Use of |linspacen| for vectors
% For (column) vector input |linspacen| produces a matrix

A=zeros(1,10)
B=ones(size(A));
C=linspacen(A(:),B(:),n)

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);
title('A surface gradient from vector A to vector B','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
surf(C);
colormap(cMap); caxis([min(C(:)) max(C(:))]);  colorbar; 
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% Use of |linspacen| for 2D matrices
% For qxr input matrices |linspacen| produces a qxrxn output matrix whereby
% the entries go from input A to B allong the last dimension in n steps

A=zeros(5,5)
B=ones(size(A));
C=linspacen(A,B,n)

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);
title('A 3D gradient from A to B','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

[Fm,Vm,Cm]=ind2patch(1:numel(C),C,'vu');

patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',1);

axis equal; view(3); axis tight; axis vis3d; grid off;
colormap(cMap); caxis([min(C(:)) max(C(:))]);  colorbar; 
set(gca,'FontSize',fontSize);
drawnow;

%% Use of |linspacen| for higher order matrices
% For qxrx... input matrices |linspacen| produces a qxrxn output matrix
% whereby the entries go from input A to B allong the last dimension in n
% steps. 

A=zeros(5,5,5);
B=ones(size(A));
C=linspacen(A,B,n);

%%
% Plotting results
hf1=figuremax(figColor,figColorDef);
subplot(1,3,1);
title('4D out put: First set','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

c=C(:,:,:,1);
[Fm,Vm,Cm]=ind2patch(1:numel(c),c,'vu');

patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',1);

axis equal; view(3); axis tight; axis vis3d; grid off;
colormap(cMap); caxis([min(C(:)) max(C(:))]); colorbar; 
set(gca,'FontSize',fontSize);

subplot(1,3,2);
title('4D out put: Intermediate set','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

c=C(:,:,:,round(size(C,4)/2));
[Fm,Vm,Cm]=ind2patch(1:numel(c),c,'vu');

patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',1);

axis equal; view(3); axis tight; axis vis3d; grid off;
colormap(cMap); caxis([min(C(:)) max(C(:))]); colorbar; 
set(gca,'FontSize',fontSize);

subplot(1,3,3);
title('4D out put: Last set','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

c=C(:,:,:,end);
[Fm,Vm,Cm]=ind2patch(1:numel(c),c,'vu');

patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',1);

axis equal; view(3); axis tight; axis vis3d; grid off;
colormap(cMap); caxis([min(C(:)) max(C(:))]); colorbar; 
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