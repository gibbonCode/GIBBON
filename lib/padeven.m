function [M_pad,L_even_dim]=padeven(M)

L_even_dim=iseven(size(M));
M_pad=zeros(size(M)+L_even_dim);
M_pad(1:1:size(M,1),1:1:size(M,2),1:1:size(M,3))=M;

