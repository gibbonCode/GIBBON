function [M_new]=imageResample(M,voxelSize,voxelSizeNew)

voxelSize=voxelSize(:)';
voxelSizeNew=voxelSizeNew(:)';

siz=size(M);
FOV=siz.*voxelSize; 
sizNew=ceil(FOV./voxelSizeNew);

[Jn,In,Kn]=meshgrid(1:1:sizNew(2),1:1:sizNew(1),1:1:sizNew(3));
[X,Y,Z]=im2cart(In(:),Jn(:),Kn(:),voxelSizeNew);
[I,J,K]=cart2im(X,Y,Z,voxelSize);
I=round(I);
J=round(J);
K=round(K);

logicValid=I<=siz(1) & J<=siz(2) & K<=siz(3);

I=I(logicValid);
J=J(logicValid);
K=K(logicValid);

IND=sub2ind(siz,I,J,K); 

M_new=nan(sizNew); 
M_new(logicValid)=M(IND);
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
