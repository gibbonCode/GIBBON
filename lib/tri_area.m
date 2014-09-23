function A=tri_area(F,V)

if size(V,2)==2
    V(:,3)=0;
end

V12=V(F(:,2),:)-V(F(:,1),:);
V13=V(F(:,3),:)-V(F(:,1),:);
A=0.5.*sqrt(sum(cross(V12,V13,2).^2,2));
