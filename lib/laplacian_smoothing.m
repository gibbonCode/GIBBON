function V=laplacian_smoothing(V,IND_V,L,n)

[I,J,v] = find(IND_V);

for i=1:n;
    Xp=accumarray({I,J},V(v,1),size(IND_V),[],NaN);
    Yp=accumarray({I,J},V(v,2),size(IND_V),[],NaN);
    Zp=accumarray({I,J},V(v,3),size(IND_V),[],NaN);
    Vp=[nanmean(Xp,2) nanmean(Yp,2) nanmean(Zp,2)];
    V=V+L.*(Vp-V);
end

end