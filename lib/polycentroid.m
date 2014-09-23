function [Xc,Yc]=polycentroid(X,Y)

%N.B. 
% Assumes row vectors or matrices whereby each row describes a polygon with
% points appearing in the order defining the polygon

meanX=mean(X,2);
meanY=mean(Y,2);
X=X-meanX*ones(1,size(X,2));
Y=Y-meanY*ones(1,size(Y,2));

A = polyarea(X',Y');
Xc=sum((X(:,1:end-1)+X(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));
Yc=sum((Y(:,1:end-1)+Y(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));

Xc=Xc+meanX;
Yc=Yc+meanY;

end