function [F_uni,V_uni,C_uni,IND_V,IND_F,F_count]=unique_patch(F,V,C,epsOrderDiff)

numFacesIni=size(F,1);

%Removing unused vertices
[F,V]=removeNotIndexed(F,V);

%Removing double vertices
[V_uni,IND_V,IND_IND]=uniqueEps(V,'rows',epsOrderDiff);

% [~,IND_V,IND_IND]=unique(pround(V,epsOrderDiff),'rows');
% V_uni=V(IND_V,:);

F=IND_IND(F); %Fix indices in F

%Removing double FACES
[F_uni,IND_F,IND_F_2]=uniqueIntegerRow(F);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Fixing face colors, shared faces now obtain mean colour
if ~isempty(C)
    sharedColourMatrixSparse=sparse(IND_F_2,1:numFacesIni,C,numFacesUni,numFacesIni,numFacesIni);    
    C_uni=full(sum(sharedColourMatrixSparse,2))./F_count;
else
    C_uni=[];
end
