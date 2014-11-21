%% box_indices
% Below is a demonstration of the features of the |box_indices| function
%

%% Syntax
% |[IND]=box_indices(siz);|

%% Description 
% The |box_indices| function returns the indices of the outer boundary
% elements of an array (which can be thought of as defining a box), i.e.
% the indices of the first and last row, columsn, slice, etc.. 

%% Examples

%%
clear; close all; clc; 

%% 
% Plot settings
fig_color='w'; fig_colordef='white'; 
faceAlpha1=1;
faceAlpha2=0.65;
edgeColor1='none';
edgeColor2='none';

%% Example: |box_indices| for 2D arrays

siz=[25 25];
M=ones(siz);
[indBox]=box_indices(size(M));
M(indBox)=0; %setting edge indices to 0 for visualization

%%
% Plotting results
figuremax(fig_color,fig_colordef);
title('Black pixels denote edge entries'); hold on; 
imagesc(M);
axis equal; axis tight; axis vis3d; grid off;  
colormap gray; caxis([0 1]); colorbar;
drawnow;

%% Example: |box_indices| for 3D arrays

siz=[25 25 25];
M=ones(siz);
[indBox]=box_indices(size(M));
M(indBox)=0; %setting edge indices to 0 for visualization

%%
% Plotting results

% Creating patch data for voxel display
logicPlot=false(size(M));
logicPlot(:,:,round(size(M,3)/2))=1; 
logicPlot(:,round(size(M,2)/2),:)=1; 
logicPlot(round(size(M,1)/2),:,:)=1; 
[F,V,C]=ind2patch(logicPlot,M,'v'); 

figuremax(fig_color,fig_colordef);
title('Black voxels denote edge entries'); hold on; 
xlabel('J - columns');ylabel('I - rows'); zlabel('K - slices'); hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','g','FaceAlpha',1);

axis equal; view(3); axis tight; axis vis3d; grid off;  
colormap gray; caxis([0 1]); colorbar;
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>