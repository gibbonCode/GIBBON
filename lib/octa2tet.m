function [Et,Vt,Ct]=octa2tet(E,V)


%%

%Get octahedron element faces
C=(1:1:size(E,1))';
[F,CF]=element2patch(E,C,'octa6');

%Get centre coordinates
Vn=patchCentre(E,V);

%Append new centre coordinates
Vt=[V;Vn]; 

%Create new tetrahedral element array 
indAdd=((size(V,1)+1:size(Vt,1))'*ones(1,8));
Et=[indAdd(:) F];

%Gather color information
Ct=CF;

