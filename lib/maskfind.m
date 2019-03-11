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
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
