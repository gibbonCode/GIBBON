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

 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
