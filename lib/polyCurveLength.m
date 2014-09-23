function D = polyCurveLength(V)

D=zeros(size(V,1),1);
D(2:end)=cumsum(sqrt(sum(diff(V,1,1).^2,2)));

