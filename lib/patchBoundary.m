function [Eb]=patchBoundary(F,V)

%Get non-unique edges
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

Es=sort(E,2);

sizVirt=size(V,1)*ones(1,2);
ind_E=sub2ind(sizVirt,Es(:,1),Es(:,2));

[~,indUni1,~]=unique(Es,'rows'); %Get indices for unique edges

E_uni=E(indUni1,:);

ind_E_uni=ind_E(indUni1,:);

ind_EF=(1:size(F,1))';
ind_EF=ind_EF(:,ones(1,size(F,2)))'; 
ind_EF=ind_EF(:);

ind=1:1:size(E,1);
ind=ind(~ismember(ind,indUni1));
ind_E_cut=ind_E(ind,:);
L=~ismember(ind_E_uni,ind_E_cut);

Eb=E_uni(L,:);
