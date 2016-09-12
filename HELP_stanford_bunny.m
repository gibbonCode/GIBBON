%% stanford_bunny
% Below is a demonstration of the |stanford_bunny| function

%%
close all; clc; clear;

% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=1;

%% THE STANFORD BUNNY
% The stanford_bunny function generates patch data (faces and vertices)
% defining a relatively coarse representation of the "Stanford bunny" which
% is a commonly used test model in computer graphics. 
% The surface is not entirely closed as can be seen in the second figure
% below. 
%
% This MATLAB implementation is based on the coarse representation
% downloadable from: 
% http://www.cc.gatech.edu/projects/large_models/bunny.html
% 
% See also:
% http://www.gvu.gatech.edu/people/faculty/greg.turk/bunny/bunny.html 
% http://graphics.stanford.edu/data/3Dscanrep/
%
% Turk G, Levoy M. Zippered polygon meshes from range images.
% Proceedings of the 21st annual conference on Computer graphics and
% interactive techniques - SIGGRAPH  �94 [Internet]. New York, New York,
% USA: ACM Press; 1994;311�8. Available from:
% http://portal.acm.org/citation.cfm?doid=192161.192241

%Obtaining patch data
[F,V]=stanford_bunny;

%Visualisation

C=rand(size(F,1),1); %random face color values

hf=figuremax(fig_color,fig_colordef); 
title('The Stanford bunny','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);

set(hp,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
colormap autumn; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting phong;

hf=figuremax(fig_color,fig_colordef); 
title('The Stanford bunny','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);

set(hp,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
colormap autumn; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  axis vis3d; axis off;
camlight('headlight'); lighting phong;
view(24.5,-44);

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
