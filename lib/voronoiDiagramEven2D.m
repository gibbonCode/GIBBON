function [V,F]=voronoiDiagramEven2D(varargin)

%% DERIVE STANDARD 2D VORONOI DIAGRAM

%Get Delaunay tesselation
DT=varargin{1};

%Computer standard Voronoi output containing vertices matrix and faces cell
[Vv,Fv_cell] =voronoiDiagram(DT);
numCell=numel(Fv_cell);

%Replace inf entries by NaN
Vv(isinf(Vv))=NaN;

%% Creating even Voronoi cell sampling

%Get number of points 
switch nargin    
    case 2
        np=varargin{2};
        interpMethod='linear';
    case 3
        np=varargin{2};
        interpMethod=varargin{3};
    otherwise
        np=[];
        interpMethod='linear';
end

%Use maximum occuring number of points if np=[]
if isempty(np)
    np=max(cellfun(@numel,Fv_cell));
end

%Resampling Voronoi cells
V=zeros(numCell*np,2); %Vertex matrix
F=reshape([1:size(V,1)]',np,size(V,1)/np)'; %Face matrix
for q=1:1:numCell
    
    %Get current Voronoi cell vertices   
    vv=Vv(Fv_cell{q},:);
    
    %Resample cell all vertices are valid
    if any(isnan(vv))
        vv=NaN;
    else        
        [vv]=evenlySampleCurve(vv,np,interpMethod,1);
    end
    
    %Store new points in array
    V(F(q,:),:)=vv;        
    
end

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
