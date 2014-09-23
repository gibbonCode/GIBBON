function [indFree]=freeBoundaryPatch(F)

%Removing double FACES
[F_uni,IND_F,IND_F_2]=uniqueIntegerRow(F);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:size(F,1),1,numFacesUni,size(F,1),size(F,1));
F_count=full(sum(logicColourMatrixEntry,2));
logicFree=F_count==1; 
indFree=IND_F(logicFree);
