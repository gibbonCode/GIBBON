function [F,V]=triMeshEquilateral(minV,maxV,pointSpacing)

% function [F,V]=triMeshEquilateral(minV,maxV,pointSpacing)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/28/04
% ------------------------------------------------------------------------
%%

%Region "Field Of View" size
FOV=abs(maxV-minV);

%Calculate number of points in x direction
nx=ceil(FOV(1)./pointSpacing)+1;
xRange=linspace(minV(1)-pointSpacing,maxV(1)+pointSpacing,nx+2);

%Actual point spacing used
pointSpacingReal_X=abs(xRange(1)-xRange(2));

%Calculate number of points in x direction
pointSpacingReal_Y=pointSpacingReal_X.*0.5*sqrt(3);
yRange=minV(2)-pointSpacingReal_Y:pointSpacingReal_Y:maxV(2)+pointSpacingReal_Y;
% yRange=[yRange(1)-pointSpacingReal_Y  yRange  yRange(end)+pointSpacingReal_Y];

%Create normal mesh
[X,Y]=meshgrid(xRange,yRange);

%Offset mesh in X direction to obtain aproximate equilateral triangular mesh
X(1:2:end,:)=X(1:2:end,:)+(pointSpacing/2);

V=[X(:) Y(:)];
V(:,3)=0; 

%Compute connectivity
[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1); %Point indices

%"Forward slash"
TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
F1 = sub2ind(size(X),TRI_I,TRI_J);

%"Backward slash"
TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
F2 = sub2ind(size(X),TRI_I,TRI_J);

F_I=[I(:);I(:)];
logicSlashType=iseven(F_I);

F=[F1(~logicSlashType,:);F2(logicSlashType,:)];

 
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
