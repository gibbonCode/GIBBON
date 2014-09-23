function [F,V]=discQuadMesh(nElements,r,f)


%% Creating central regular quad mesh

nElements=nElements+~iseven(nElements);%Force even

[X_centralMesh,Y_centralMesh]=meshgrid(linspace(-1,1,nElements+1));
[F_centralMesh,V_centralMesh] = surf2patch(X_centralMesh,Y_centralMesh,zeros(size(X_centralMesh)));
V_centralMesh=V_centralMesh(:,1:2);

%Edge of central mesh
logicCentralMeshEdge=(X_centralMesh==1)|(Y_centralMesh==1)|(X_centralMesh==-1)|(Y_centralMesh==-1);
nEdge=(nElements*4);

% Scaling radius
[ThetaMesh,RadiusMesh]=cart2pol(V_centralMesh(:,1),V_centralMesh(:,2));
RadiusMesh=f*(1/2)*sqrt(2)*RadiusMesh;
[V_centralMesh(:,1),V_centralMesh(:,2)]=pol2cart(ThetaMesh,RadiusMesh);

%% Creating outer mesh

RadiusOuterEdge=ones(1,nEdge);
ThetaOuterEdge=linspace(0,pi*2,nEdge+1); 
ThetaOuterEdge=ThetaOuterEdge(2:end)-pi;

[xOuterEdge,yOuterEdge]=pol2cart(ThetaOuterEdge,RadiusOuterEdge);
V_outerEdge=[xOuterEdge(:) yOuterEdge(:)];

V_innerEdge=V_centralMesh(logicCentralMeshEdge,:);
[ThetaEdge,RadiusEdge]=cart2pol(V_innerEdge(:,1),V_innerEdge(:,2));
[ThetaEdge,sortInd]=sort(ThetaEdge);
RadiusEdge=RadiusEdge(sortInd);
[V_innerEdge(:,1),V_innerEdge(:,2)]=pol2cart(ThetaEdge,RadiusEdge);

[Xr]=linspacen(V_innerEdge(:,1),V_outerEdge(:,1),nElements/2+1); Xr(end+1,:)=Xr(1,:);
[Yr]=linspacen(V_innerEdge(:,2),V_outerEdge(:,2),nElements/2+1); Yr(end+1,:)=Yr(1,:);

logicEdge=false(size(Xr));logicEdge(:,end)=1;
[Fs2,Vs2,logicEdge] = surf2patch(Xr,Yr,zeros(size(Xr)),logicEdge);
Vs2=Vs2(:,1:2);

V=[V_centralMesh;Vs2];
F=[F_centralMesh;Fs2+size(V_centralMesh,1)];
logicEdge=[false(size(V_centralMesh,1),1); logicEdge];

%% Removing double points

[~,IND_V,IND_IND]=unique(pround(V,5),'rows'); %N.B. note 5th decimal rounding used here
V=V(IND_V,:);
F=IND_IND(F);

%Scaling radius
[ThetaMesh,RadiusMesh]=cart2pol(V(:,1),V(:,2));
RadiusMesh=r*RadiusMesh;
[V(:,1),V(:,2)]=pol2cart(ThetaMesh,RadiusMesh);

end



