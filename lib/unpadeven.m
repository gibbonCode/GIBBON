function [M]=unpadeven(M_pad,L_even_dim)

num_dims=ndims(M_pad);
siz=size(M_pad)-L_even_dim;
switch num_dims
    case 2
        M=M_pad(1:1:siz(1),1:1:siz(2));
    case 3
        M=M_pad(1:1:siz(1),1:1:siz(2),1:1:siz(3));        
end




