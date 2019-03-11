function [F,V]=patchNanFix(F,V)

numPoints=size(V,1); %Get original number of entries in vertex list

logicValid=~any(isnan(V),2); %Logic for valid vertices
logicValid_F=all(logicValid(F),2); %Logic vor valid faces

F=F(logicValid_F,:); %Get valid faces
V=V(logicValid,:); %Get valid vertices

%Fix indices in faces matrix
indFix1=1:nnz(logicValid);
indFix2=zeros(numPoints,1);
indFix2(logicValid)=indFix1;
F=indFix2(F); %Fix indices
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
