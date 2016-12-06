function ind=subMat2ind(siz,I)

k=cumprod(siz);
ind = I(:,1);
for q=2:1:size(I,2)
    ind = ind + (I(:,q) - 1).*k(q-1);
end
