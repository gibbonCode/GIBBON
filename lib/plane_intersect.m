function X=plane_intersect(V1,V2,V3,N1,N2,N3)

%%

DET=ndet(N1,N2,N3);
DET(DET==0)=NaN;

%  X= (1./DET).*[(dot(V(:,1),N(:,1)).*cross(N(:,2),N(:,3))) + (dot(V(:,2),N(:,2)).*cross(N(:,3),N(:,1))) + (dot(V(:,3),N(:,3)).*cross(N(:,1),N(:,2)))]
 
X= ((1./DET)*ones(1,3)).* (...
     ((dot(V1,N1,2)*ones(1,3)).*cross(N2,N3,2)) + ...
     ((dot(V2,N2,2)*ones(1,3)).*cross(N3,N1,2)) + ...
     ((dot(V3,N3,2)*ones(1,3)).*cross(N1,N2,2))...
     );


end
 
 