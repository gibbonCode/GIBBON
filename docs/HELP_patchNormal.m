%% patchNormal
% Below is a demonstration of the features of the |patchNormal| function

%%
clear; close all; clc;

%% Syntax
% |[N,Vn,Nv]=patchNormal(F,V);|

%% Description 
% Compute the surface normals N for the input path data defined by the
% faces F and vertices V. Optional additional outputs include the face
% centres Vn or the vertex normals Nv. 

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

[N,P,NV]=patchNormal(F,V);

%%
% Visualization

s=mean(patchEdgeLengths(F,V));

cFigure;
subplot(1,2,1);
title('Face normals')
gpatch(F,V,'bw');
quiverVec(P,N,s,'k')
axisGeom; camlight headlight;

subplot(1,2,2);
title('Vertex normals')
gpatch(F,V,'bw');
quiverVec(V,NV,s,'k')
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
