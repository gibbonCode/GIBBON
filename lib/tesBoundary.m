function [indBounary]=tesBoundary(F,V)

if numel(V)==1
    numPoints=V;
else
    numPoints=size(V,1);
end

Fbs=sort(F,2);

sizVirt=numPoints*ones(1,size(F,2));

ind_F=subMat2ind(sizVirt,Fbs);

[~,indUni1,~]=unique(Fbs,'rows'); %Get indices for unique faces
F_uni=F(indUni1,:);

ind_F_uni=ind_F(indUni1,:);

ind=1:1:size(F,1);
ind=ind(~ismember(ind,indUni1));
ind_Fb_cut=ind_F(ind,:);
L_uni=~ismember(ind_F_uni,ind_Fb_cut);

indBounary=indUni1(L_uni,:);
