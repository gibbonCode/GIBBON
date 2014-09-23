function [I,J,K,IND]=logic2subind(L)

[J,I,K]=meshgrid(1:1:size(L,2),1:1:size(L,1),1:1:size(L,3));
I=I(L); J=J(L); K=K(L); 
IND = sub2ind(size(L),I,J,K);
