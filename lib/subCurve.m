function [Vs]=subCurve(V,n)

if n==0 %No subdevision
    Vs=V;     
elseif n>1 %Subdevision of segments
    Vs=zeros(size(V,1)+(size(V,1)-1)*n,size(V,2));
    for q=1:1:size(V,2);
        X=V(:,q);
        XX=linspacen(X(1:end-1),X(2:end),n+2);
        XX=XX(:,1:end-1)';
        Vs(1:end-1,q)=XX(:);
        Vs(end,q)=V(end,q);
    end
elseif n<1 %Throw error
    error('n should be >=0')
end
