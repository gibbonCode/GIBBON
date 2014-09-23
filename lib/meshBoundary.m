function [Fb,F_count]=meshBoundary(E)

%Get mesh faces
[F,~]=element2patch(E,[]);

%Check shared faces faces
numFacesIni=size(F,1);
[F_uni,~,IND_F_2]=uniqueIntegerRow(F);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Compose boundary set from faces that are used once
Fb=F_uni(F_count==1,:);
