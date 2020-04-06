function D=edgeLengths(E,V)

D=sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));