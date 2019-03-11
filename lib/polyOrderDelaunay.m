function [V_fix,indFix]=polyOrderDelaunay(V,fd)

%%
% N.B. FUNCTION UNFINISHED AND NOT FULLY VALIDATED

%%
% Get Delaunay tesselation of points set
DT=delaunayTriangulation(V);
% Vd=DT.Points;
nPointsOriginal=size(V,1);

% Add Voronoi vertices to set
Vv = DT.voronoiDiagram();
Vv=Vv(any(~isinf(Vv),2),:); %Remove infinity vertex
DT.Points=[DT.Points;Vv];

% The Delaunay edges that connect pairs of sample points represent the
% boundary.
edgeInd = DT.edges();
logicBoundaryEdges=all(edgeInd<=nPointsOriginal,2);
E = edgeInd(logicBoundaryEdges,:);

%Removing invalid edges
D=sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
L_invalid=(D./mean(D(:)))>fd;
E=E(~L_invalid,:);

indUni=unique(E(:));
indFix=nan(numel(indUni),1);

Ec=E;
sizE=size(E);

%% Loop to get path

q=1;
indNow=indUni(1);
% hp=[];
while 1    
    if q==numel(indUni)+1
        break
    end
    
    indFix(q)=indNow;
    indE=find(Ec==indNow,1);
    [rowE,~] = ind2sub(sizE,indE);
    Ec(rowE,:)=NaN;
    indV=E(rowE,:);    
    indNow=indV(~ismember(indV,indFix));        
 
    q=q+1;        
%     V_fix=V(indFix(~isnan(indFix)),:);
%     delete(hp);    
%     hp=plotV(V_fix, 'k+-','MarkerSize',markerSize1);
%     drawnow; pause(0.01);
end
V_fix=V(indFix,:);

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
