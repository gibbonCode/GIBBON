function [Es,Vs]=subTet(E,V,splitMethod)

% function [Es,Vs]=subTet(E,V,splitMethod)
% ------------------------------------------------------------------------
% Sub devides the input tetrahedral mesh into a denser mesh by splitting
% the elements. Two methods are available. If splitMethod==1 then the faces
% are split similar to what is expected for [Fs,Vs]=subtri(F,V,1,1). I.e.
% for side faces mid-edge nodes are introduced so edges are split. This 
% yields 4 new corner tetrahedral elements which could inherit shape
% quality from the original tetrahedral element. However the entire
% tetrahedron cannot me split this way and a central octahedron remains.
% This octahedron is split using its Delaunay tesselation representation
% and yields 4 new tetrahedrons which differ in shape quality from the
% original tetrahedron. If splitMethod==2 the tetrahedron is split into 4
% tetrahedra by introducing a single new central node for each element and
% by connecting the side faces
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/03/04 
% ------------------------------------------------------------------------

switch splitMethod
    case 1 %Split faces + central delaunay of octahedron
        edgeMat=[E(:,[1 2]); E(:,[2 3]);  E(:,[3 1]); E(:,[1 4]); E(:,[3 4]); E(:,[2 4])]; %Edges matrix
        E_sort=sort(edgeMat,2); %Sorted edges matrix
        [~,ind1,~]=unique(E_sort,'rows');
        edgeMat=edgeMat(ind1,:);
        
        numPoints = size(V,1);
        numEdges = size(edgeMat,1);
        
        % Get indices of the three edges associated with each face
        A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
        A = max(A,A'); %Copy symmetric
        
        %Indices for A matrix
        indA_12=E(:,1)+(E(:,2)-1)*numPoints;
        indA_23=E(:,2)+(E(:,3)-1)*numPoints;
        indA_31=E(:,3)+(E(:,1)-1)*numPoints;
        indA_14=E(:,1)+(E(:,4)-1)*numPoints;
        indA_24=E(:,2)+(E(:,4)-1)*numPoints;
        indA_34=E(:,3)+(E(:,4)-1)*numPoints;
        
        %Get indices for vertex array
        indV_12=full(A(indA_12));
        indV_23=full(A(indA_23));
        indV_31=full(A(indA_31));
        indV_14=full(A(indA_14));
        indV_24=full(A(indA_24));
        indV_34=full(A(indA_34));
        
        %Create element array
        Es=[E(:,1)  indV_12 indV_31 indV_14;... %Corner tet 1
            indV_12 indV_23 indV_24 E(:,2);... %Corner tet 2
            E(:,3)  indV_31 indV_23 indV_34;... %Corner tet 3
            indV_14 indV_24 indV_34 E(:,4);... %Corner tet 4
            indV_34 indV_31 indV_23 indV_24;... %Octahedral tet 1
            indV_31 indV_12 indV_23 indV_24;... %Octahedral tet 2
            indV_31 indV_12 indV_24 indV_14;... %Octahedral tet 3
            indV_34 indV_31 indV_24 indV_14;... %Octahedral tet 4
            ];
        
        %Create vertex array
        Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
        Vs = [V; Vn]; %Join point sets
        
            
    case 2 %Add central point and connect all faces to central point to yield elements
        %Get faces
        [F,faceInd]=element2patch(E,(1:size(E,1))');    
        
        %Create element array
        Es=[F faceInd+size(V,1)];
        
        %Create vertex array
        X=V(:,1); Y=V(:,2); Z=V(:,3);
        if size(E,1)==1 %Indexing behaviour (e.g. X(E)) differs for a single element 
            Vn=[mean(X(E),1) mean(Y(E),1) mean(Z(E),1)]; %new mid-element points
        else
            Vn=[mean(X(E),2) mean(Y(E),2) mean(Z(E),2)]; %new mid-element points
        end
        Vs = [V; Vn]; %Join point sets
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
