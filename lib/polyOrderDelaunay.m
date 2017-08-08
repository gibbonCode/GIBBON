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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
