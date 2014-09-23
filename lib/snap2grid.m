function [Xs,Ys,Zs]=snap2grid(X,Y,Z,c,v)

Xs=c(1)+round((X-c(1))./v(1)).*v(1);
Ys=c(2)+round((Y-c(2))./v(2)).*v(2);
Zs=c(3)+round((Z-c(3))./v(3)).*v(3);

end



