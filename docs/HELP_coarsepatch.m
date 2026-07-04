%% coarsepatch
% Below is a demonstration of the features of the |coarsepatch| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,Eb]=coarsepatch(F,V,f);|

%% Description 
% This function can reduce the density of patch data defined by faces |F|
% and vertices |V|. The input |f| either defines a scaling factor between
% 0-1 or a desired number of faces. 
% This function uses MATLAB's |[F,V] = reducepatch(F,V,f,'fast');|. However
% it also proceeds to clean up the mesh to avoid several reducepatch issues
% (collapsed triangles etc). 
% The 3rd optional output consists of the boundary edges |Eb|. If the input
% does not have boundary edges then the coarsening operation introduced
% holes. 

%% Examples 
% 

%% Specifying a scaling factor
[F1,V1]=graphicsModels(9);

f=1/3; %Scaling factor (if <=1) or new number of faces

[F2,V2]=coarsepatch(F1,V1,f);

%%

cFigure; 
subplot(1,2,1);
title(['Original: ',num2str(size(V1,1)),' vertices, ',num2str(size(F1,1)),' faces']);
gpatch(F1,V1,'w','k');
axisGeom;
camlight headlight;
view(0,0);
axis off; 

subplot(1,2,2);
title(['Resampled: ',num2str(size(V2,1)),' vertices, ',num2str(size(F2,1)),' faces']);
gpatch(F2,V2,'w','k');
axisGeom;
camlight headlight;
view(0,0);
axis off; 

drawnow; 

%% Specifying desired number of faces

[F1,V1]=graphicsModels(9);

f=1500; %Scaling factor (if <=1) or new number of faces

[F2,V2]=coarsepatch(F1,V1,f);

%%

cFigure; 
subplot(1,2,1);
title(['Original: ',num2str(size(V1,1)),' vertices, ',num2str(size(F1,1)),' faces']);
gpatch(F1,V1,'w','k');
axisGeom;
camlight headlight;
view(0,0);
axis off; 

subplot(1,2,2);
title(['Resampled: ',num2str(size(V2,1)),' vertices, ',num2str(size(F2,1)),' faces']);
gpatch(F2,V2,'w','k');
axisGeom;
camlight headlight;
view(0,0);
axis off; 

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
