function [h]=surf_slice(M,X,Y,Z,S)

h=[ ];
for i=S
    hi=surf(X(:,:,i), Y(:,:,i),Z(:,:,i),M(:,:,i),'EdgeColor','none'); 
    h=[h hi];
end
