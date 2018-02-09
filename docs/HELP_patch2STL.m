%% patch2STL
% Below is a demonstration of the features of the |patch2STL| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white'; 
edgeColor1='none';
edgeColor2='none';

%% 
% Get patch data
[F,V]=stanford_bunny;

%%
% Plotting the model 

figuremax(fig_color,fig_colordef);
xlabel('X');ylabel('Y'); zlabel('Z'); hold on;

patch('Faces',F,'Vertices',V,'FaceColor','b','EdgeColor','k','FaceAlpha',1);

view(3); axis equal; axis tight; axis vis3d; grid on; 
camlight('headlight');
lighting phong; axis off; 
drawnow;

%% Exporting an STL file from patch data

%Set main folder and fileName
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'standford_bunny.stl'); 

patch2STL(fileName,V,F,[],'standford_bunny');

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
