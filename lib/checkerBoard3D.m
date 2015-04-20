function M=checkerBoard3D(siz)

[I,J,K]=ndgrid(1:1:siz(1),1:1:siz(2),1:1:siz(3));
logic_ij=((iseven(I)| iseven(J)) & ((iseven(I)~=iseven(J))));
M=false(siz);
M(iseven(K))=logic_ij(iseven(K));
M(~iseven(K))=~logic_ij(~iseven(K));