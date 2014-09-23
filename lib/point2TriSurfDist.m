function D=point2TriSurfDist(F1,V1,V2)

%Get face normals
[N1,~]=trinorm(F1,V1);

%Compute face centres
X=V1(:,1); Y=V1(:,2); Z=V1(:,3);
XF=X(F1); YF=Y(F1); ZF=Z(F1);

%Position vectors for face centres
v1=[mean(XF,2) mean(YF,2) mean(ZF,2)];

%Find closest triangles
Dm=dist(V2,v1');
[~,indMin]=min(Dm,[],2); %Get indices of closest

%Difference vectors between face centres and points
x2=V2-v1(indMin,:);

%Find distance to closest triangles
D=abs(dot(N1(indMin,:),x2,2));