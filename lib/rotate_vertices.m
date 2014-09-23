function [Vr]=rotate_vertices(V,R_mat)

Vr=zeros([size(V) size(R_mat,3);]);

for i=1:1:size(R_mat,3);
    
    R=R_mat(:,:,i);
    Vr(:,:,i)=(R*V')';
    
end

end