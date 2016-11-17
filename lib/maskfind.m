function [IND_mask]=maskfind(M,IND,MASK_I,MASK_J,MASK_K)

% function [IND_mask]=maskfind(M,IND,I,J,K)
% ------------------------------------------------------------------------
% This function finds the indices of the elements found inside the mask
% defined by MASK_I, MASK_J and MASK_K. The matrix IND_mask is numel(IND) x
% numel(MASK_I) in size. N.B. zeros are used when indices are not found. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 07/04/2010
% ------------------------------------------------------------------------

%%

%Get image size
siz=size(M);
if size(M,3)==1
    siz(3)=1;
end

%Force row vector input
MASK_I=MASK_I(:)';
MASK_J=MASK_J(:)';
MASK_K=MASK_K(:)';

%Get subscript indices for IND
[Im,Jm,Km]=ind2sub(siz,IND(:));

%Create adjacency subscript indices
Im=Im(:,ones(1,numel(MASK_I)))+MASK_I(ones(numel(Im),1),:);
Jm=Jm(:,ones(1,numel(MASK_J)))+MASK_J(ones(numel(Jm),1),:);
Km=Km(:,ones(1,numel(MASK_K)))+MASK_K(ones(numel(Km),1),:);

%Get logic for valid indices. 
L_valid= ~((Im<1 | Im>siz(1)) | (Jm<1 | Jm>siz(2)) | (Km<1 | Km>siz(3)));

%Get linear indiced of adjacenty indices
[ind_m]=sub2ind(siz,Im(L_valid),Jm(L_valid),Km(L_valid));

IND_mask=zeros(size(Im));
IND_mask(L_valid)=ind_m;

