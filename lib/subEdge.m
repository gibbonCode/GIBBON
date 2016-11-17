function [F_add,V_add]=subEdge(F,V,n)

%Edges matrix
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

V_add=zeros(size(E,1)*(n+1),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    XE=X(E);
    X_add=linspacen(XE(:,1),XE(:,2),n+2);  
    X_add=X_add(:,1:end-1)';
    V_add(:,q)=X_add(:);
end

ind=(1:1:size(V_add,1))';
F_add=reshape(ind,[(n+1)*size(F,2),size(V_add,1)./((n+1)*size(F,2))])';

%Merge non-unique nodes
numDigitKeep=5; 
% [~,ind1,ind2]=unique(pround(V_add,numDigitKeep),'rows');
try
    [~,ind1,ind2]=unique(round(V_add,numDigitKeep,'significant'),'rows');
catch
    [~,ind1,ind2]=unique(sround(V_add,numDigitKeep),'rows');
end

V_add=V_add(ind1,:);
F_add=ind2(F_add);