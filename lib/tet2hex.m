function [HEX,Vhex]=tet2hex(TET,V)

%% Get faces
F=[TET(:,[2 1 3]); TET(:,[1 2 4]); TET(:,[2 3 4]); TET(:,[3 1 4])];

%% Deriving new coordinate sets

%The original vertices
X=V(:,1); Y=V(:,2); Z=V(:,3);

%The mid-face points
Vmf=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];

%The mid points
if size(TET,1)==1
    Vm=[mean(X(TET)) mean(Y(TET)) mean(Z(TET))];
else
    Vm=[mean(X(TET),2) mean(Y(TET),2) mean(Z(TET),2)];
end

%The mid edge points
indEdge1=[1 2]; indEdge2=[2 3]; indEdge3=[3 1];
indEdge4=[1 4]; indEdge5=[2 4]; indEdge6=[3 4];
edgeIndAll=[TET(:,indEdge1);TET(:,indEdge2);TET(:,indEdge3);...
    TET(:,indEdge4);TET(:,indEdge5);TET(:,indEdge6);];

if size(TET,1)==1
    Xt=reshape(X(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
    Yt=reshape(Y(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
    Zt=reshape(Z(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
else
    Xt=X(edgeIndAll);
    Yt=Y(edgeIndAll);
    Zt=Z(edgeIndAll);
end
Vme=[mean(Xt,2) mean(Yt,2) mean(Zt,2)];

%the joint coordinate set
Vhex=[V; Vme; Vmf; Vm];

%% Creating hexahedral elements

%Numbers for fixing point indices
numVert=size(V,1);
numTET=size(TET,1);
numVme=size(Vme,1);
numVmf=size(Vmf,1);

%Indices for corner, edge-, face- points etc.
mixInd=[1 1 1 3 4 2 0 4;...
        2 2 1 1 5 3 0 2;...
        3 3 1 2 6 4 0 3;...
        4 4 2 5 6 4 0 3];

HEX=zeros(size(TET,1).*4,8); %The matrix for the hexahedral element spec.
for hexInd=1:1:4;
    
    %Corner
    ind1=TET(:,mixInd(hexInd,1));
    
    %Mid edge point
    edgeInd=mixInd(hexInd,2);
    ind2=numVert+((edgeInd-1).*numTET)+(1:numTET);
    
    %Mid face point
    faceInd=mixInd(hexInd,3);
    ind3=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
    
    %Mid edge point
    edgeInd=mixInd(hexInd,4);
    ind4=numVert+((edgeInd-1).*numTET)+(1:numTET);
    
    %Mid edge point
    edgeInd=mixInd(hexInd,5);
    ind5=numVert+((edgeInd-1).*numTET)+(1:numTET);
    
    %Mid face point
    faceInd=mixInd(hexInd,6);
    ind6=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
    
    %Mid point
    ind7=numVert+numVme+numVmf+(1:numTET);
    
    %Mid face point
    faceInd=mixInd(hexInd,8);
    ind8=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
    
    h=[ind1(:) ind2(:) ind3(:) ind4(:) ind5(:) ind6(:) ind7(:) ind8(:)];       
    
    if hexInd==4
        h=h(:,[5 6 7 8 1 2 3 4]);
    end
    startInd=1+(hexInd-1).*(size(TET,1));
    endInd=startInd-1+(size(TET,1));
    HEX(startInd:endInd,:)=h;
    
end

%% Removing double VERTICES
epsOrderDiff=5;
[Vhex,~,IND_IND]=uniqueEps(Vhex,'rows',epsOrderDiff);
HEX=IND_IND(HEX);

end
