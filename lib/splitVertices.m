function [Fn,Vn]=splitVertices(F,V)

%%

%If vertices are shared they need to be disconnected by adding vertices

Vn=zeros(numel(F),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q); %The coordinate set
    XF=X(F); %The coordinates for each face vertex
    XFt=XF'; %transpose to prepare for column array
    Vn(:,q)=XFt(:); %Add column as part of new vertex matrix
end

%Fix indices in face matrix
indVn=1:numel(F);
Fn=reshape(indVn(:),size(F,2),size(F,1))';

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
