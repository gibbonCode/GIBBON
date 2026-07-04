%% patchNormPlot
% Below is a demonstration of the features of the |patchNormPlot| function

%%
clear; close all; clc;

%% Syntax
% |[hp]=patchNormPlot(F,V,a,pathType);|

%% Description 
% Visualizes surface normals for all faces F with vertices V.

%% Examples 
% 

testCase=2;
switch testCase
    case 1 %Single element square 1x1
        z=2;
        V=[0 0 0; 1 0 0; 0 1 z; 1 1 0];
        F=[1 2 4 3];            
    case 2 %Sphere triangles
        r=1;
        n=2;
        [F,V]=geoSphere(n,r);        
end

%%
% Visualization of face normals


cFigure;
gpatch(F,V,'w');

hp=patchNormPlot(F,V); %Visualize face normals

axisGeom; camlight headlight;
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
