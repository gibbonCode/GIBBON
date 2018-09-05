function [F_quad,V_quad]=tri2quadGroupSplit(F_tri,V_tri,maxAngleDeviation)

E=patchEdges(F_tri,1);
[C]=patchConnectivity(F_tri,V_tri);
EF=C.edge.face;
logicValid=all(EF>0,2);
EF=EF(logicValid,:);
E=E(logicValid,:);

FS=[F_tri(EF(:,1),:) F_tri(EF(:,2),:)];

I=(1:1:size(E,1))';
I=I(:,ones(1,6));
I=I(:);

S=sparse(I,FS,1,size(E,1),size(V_tri,1),numel(FS))>0;
[I,J]=find(S);
S=double(S);
S(S>0)=J;
S=sort(S,2,'descend');
S=full(S(:,1:4));
% L=true(size(FS));
% for q=1:1:6
%     L(:,q)=any(FS(:,q*ones(1,2))==E,2);
% end
% FS(L)=0
FQ=S;

A1=sum(patchEdgeAngles(FQ,V_tri),2);
A2=sum(patchEdgeAngles(FQ(:,[2 1 3 4]),V_tri),2);
FQ=[FQ(A1>=A2,:);FQ(A1<A2,[2 1 3 4]);];

A1=sum(patchEdgeAngles(FQ,V_tri),2);
A2=sum(patchEdgeAngles(FQ(:,[1 3 2 4]),V_tri),2);
FQ=[FQ(A1>=A2,:);FQ(A1<A2,[1 3 2 4]);];

N=patchNormal(FQ,V_tri);
FQ(N(:,3)<0,:)=fliplr(FQ(N(:,3)<0,:));

A=patchEdgeAngles(FQ,V_tri);
A(any(abs(A-(pi/2))>maxAngleDeviation,2),:)=NaN;
F_quad_sub=[];
while 1
    
    [~,indMin]=min(sum(abs(A-(pi/2)),2));
    
    A(indMin)=NaN;
    fq_now= FQ(indMin,:);
    F_quad_sub=[fq_now; F_quad_sub];
    
    logicRemove=sum(ismember(FQ,fq_now),2)>2;
    A(logicRemove,:)=NaN;
    if nnz(isnan(A))==numel(A)
        break
    end
end

F_quad_sub_tri=[F_quad_sub(:,[1 2 3]); F_quad_sub(:,[3 4 1]); F_quad_sub(:,[4 1 2]); F_quad_sub(:,[1 2 4]);];

E1=sort(patchEdges(F_tri,0),2);
indEdges_1=sub2indn(size(V_tri,1)*ones(1,2),E1);
e=reshape(indEdges_1,size(F_tri,2),size(F_tri,1))';

E2=sort(patchEdges(F_quad_sub_tri,0),2);
indEdges_2=sub2indn(size(V_tri,1)*ones(1,2),E2);
e_quad_tri=reshape(indEdges_2,size(F_quad_sub_tri,2),size(F_quad_sub_tri,1))';

logicRemove=all(ismember(e,e_quad_tri),2);
F_tri=F_tri(~logicRemove,:);

[F_quad_sub,V_quad_sub]=subQuad(F_quad_sub,V_tri,1);
[F_quad2,V_quad2]=tri2quad(F_tri,V_quad_sub);

[F_quad,V_quad]=joinElementSets({F_quad_sub,F_quad2},{V_quad_sub,V_quad2});
[F_quad,V_quad]=patchCleanUnused(F_quad,V_quad);
[F_quad,V_quad]=mergeVertices(F_quad,V_quad);
