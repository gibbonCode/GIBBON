clear; close all; clc;

%% Create example geometry

r=1;
[F,V]=hemiSphereMesh(2,r,0);
Eb=patchBoundary(F);
VF=patchCentre(F,V);
logicKeep=VF(:,3)<r/2;
F=F(logicKeep,:);
[F,V,indFix]=patchCleanUnused(F,V);
Eb=indFix(Eb);

numEdgesNeeded=size(Eb,1)+5; 

%% 
% Split edges until number of desired edges is reached

[Fn,Vn,Ebn]=triSurfSplitBoundary(F,V,Eb,numEdgesNeeded);

%%
% Visualize result

cFigure; 
subplot(1,2,1);
gpatch(F,V,'bw','k');
gpatch(Eb,V,'none','g',1,3);
axisGeom;
camlight headlight; 

subplot(1,2,2);
gpatch(Fn,Vn,'gw','k');
gpatch(Ebn,Vn,'none','b',1,3);
axisGeom;
camlight headlight; 

drawnow; 


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
