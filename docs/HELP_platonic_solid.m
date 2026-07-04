%% platonic_solid
% Below is a demonstration of the features of the |platonic_solid| function

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=2;
markerSize=5;

%%

hf=cFigure; % Open figure for plotting
pColor=gjet(5);
for q=1:1:5
    %Defining the faces (F) and vertices (V) of a platonic solid
    [V,F]=platonic_solid(q,1); %q indicates solid type, r is the radius
    
    subplot(2,3,q); hold on;
    
    hp=gpatch(F,V,pColor(q,:),'k',faceAlpha,3);
    patchNormPlot(F,V);
    
    axisGeom(gca,fontSize);    
    camlight('headlight'); lighting flat;
    axis off;
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
