function Nn=vec_normalize(N)

%Normalizing vector length
Nn=N./(sqrt(sum(N.^2,2))*ones(1,size(N,2)));

