function [F_uni,F_count]=patchCount(F)

%Removing double FACES
Fs=sort(F,2); %Sort so faces with same nodes have the same rows
[~,IND_F,IND_F_2]=unique(Fs,'rows');
F_uni=F(IND_F,:);

%Get face counts
numF=size(Fs,1);
numFuni=size(F_uni,1);
logicColourMatrixEntry=sparse(IND_F_2,1:numF,1,numFuni,numF,numF);
F_count=full(sum(logicColourMatrixEntry,2));

