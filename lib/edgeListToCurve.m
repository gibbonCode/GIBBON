function [indList]=edgeListToCurve(E)

indUni=unique(E(:));

ind1=E(1,1);
ind2=E(1,2);
indList=ones(1,numel(indUni));
indList(1)=ind1;
indList(2)=ind2;
q=3;
Es=E;
Es(1,:)=NaN;
while 1
    [indE_next,~]=find(Es==ind2,1);
    E_now=Es(indE_next,:);
    Es(indE_next,:)=NaN;
    ind3=E_now(E_now~=ind2);
    indList(q)=ind3;
    ind2=ind3;
    if nnz(isnan(Es))==numel(Es)
        break
    end
    q=q+1; 
end

end
