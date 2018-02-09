%% import_off
% Below is a demonstration of the features of the |import_off| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; 
fig_colordef='white'; 
faceAlpha=1;
fontSize=10; 

%% Import OFF file as patch data

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','OFF'); 
testCase=1; 
switch testCase
    case 1
        offName='elephant-50kv.off';
end
fileName=fullfile(pathName,offName); 

%%
% The |import_off| function imports <http://www.geomview.org> type .off
% files specifying a surface mesh (provided all faces are of the same
% type!). The output consists of (patch data) faces F and vertices V. 

[F,V] = import_off(fileName);

%%
% Plotting the models 

figuremax(fig_color,fig_colordef);
title('Imported patch data from OFF file','fontSize',fontSize);
xlabel('X','fontSize',fontSize);ylabel('Y','fontSize',fontSize); zlabel('Z','fontSize',fontSize); hold on;

patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'EdgeColor','k','FaceAlpha',faceAlpha);

view(3); axis equal; axis tight; axis vis3d; grid on; axis off; 
view(194.5,-40);
camlight('headlight');
lighting flat;
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
