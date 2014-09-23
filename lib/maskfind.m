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
siz=size(M);
if size(MASK_I,1)>1; MASK_I=MASK_I'; MASK_J=MASK_J'; MASK_K=MASK_K'; end

[Im,Jm,Km]=ind2sub(siz,IND);
Im=(Im*ones(1,numel(MASK_I)))+(ones(numel(Im),1)*MASK_I);
Jm=(Jm*ones(1,numel(MASK_J)))+(ones(numel(Jm),1)*MASK_J);
Km=(Km*ones(1,numel(MASK_K)))+(ones(numel(Km),1)*MASK_K);
L_valid= ~((Im<1 | Im>siz(1)) | (Jm<1 | Jm>siz(2)) | (Km<1 | Km>siz(3)));
[ind_m]=sub2ind(siz,Im(L_valid),Jm(L_valid),Km(L_valid));
IND_mask=zeros(size(Im));
IND_mask(L_valid)=ind_m;
