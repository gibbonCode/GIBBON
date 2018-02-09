%% graphicsModels
% Below is a demonstration of the features of the |graphicsModels| function

%% Syntax
% |[F,V]=graphicsModels(modelId);|

%% Description 
% 
%% Examples 

%%
clear; close all; clc;

% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=0.5;

%% 

hf=cFigure; 

cMap=gjet(7);
cNames={'Stanford bunny','Utah teapot','cow','parasaurolophus','femur','hip implant','elephant'};
for q=1:1:7
   subplot(3,3,q); 
    [F,V]=graphicsModels(q);
    
    title(cNames{q},'FontSize',fontSize);
    
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    
    hp=patch('Faces',F,'Vertices',V,'FaceColor',cMap(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','none');
    
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  axis vis3d; axis off;
    camlight('headlight'); lighting flat;
end

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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
