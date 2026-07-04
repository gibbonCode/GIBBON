function [varargout]=edgeVec(E,V)

% function [N,Vp,Nv]=edgeVec(E,V)
% ------------------------------------------------------------------------
% Computes the edge vectors for the input edges defined by the edge array E
% and the vertices V. The ouput is the allong edge unit vectors N. Other
% optional outputs include the edge vector origins Vp, and the vertex edge
% unit vectors. The latter are an average of both edges connected to the
% vertex. 
%
% ------------------------------------------------------------------------

%% Compute edge vectors
Vp=V(E(:,1),:); %Edge vector origin
N=V(E(:,2),:)-Vp; %Edge vector

%% Collect output
varargout{1}=N;
varargout{2}=Vp;

if nargout==3
    Nv=faceToVertexMeasure(E,V,N); %Edge vectors at vertices
    varargout{3}=Nv;
end

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
