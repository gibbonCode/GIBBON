function [VI,L]=TriScatteredInterp_nat_near(X,V,XI)

%Natural neighbour interpolation
F=scatteredInterpolant(X,V,'natural');
VI=F(XI);

%Use nearest neighbour for points outside of convex-hull
L=any(isnan(VI),2);
if any(L)
    F=scatteredInterpolant(X,V,'nearest');
    VI(L,:)=F(XI(L,:));
end











