%% stanford_bunny
% Below is a demonstration of the features of the |stanford_bunny| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=stanford_bunny(modelID)|

%% Description 
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
% interactive techniques - SIGGRAPH  -94 [Internet]. New York, New York,
% USA: ACM Press; 1994;311-8. Available from:
% http://portal.acm.org/citation.cfm?doid=192161.192241

%% Examples 
% 

%%
% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=1;

%% Obtaining the Stanford bunny patch data (default)
[F,V]=stanford_bunny;

%%
% Visualisation

cFigure;
title('The Stanford bunny','FontSize',fontSize);
gpatch(F,V,'gw');
axisGeom(gca,fontSize);
camlight('headlight');
drawnow; 

%% Obtaining the Stanford bunny patch data (other)

hf=cFigure;
gtitle('The Stanford bunny');
drawnow; 
c=gjet(6);
for q=0:5
    [F,V]=stanford_bunny(q);

    subplot(2,3,q+1);
    title(['Model ID ',num2str(q)]);
    gpatch(F,V,c(q+1,:));
    axisGeom(gca,fontSize);
    camlight('headlight');
end
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
