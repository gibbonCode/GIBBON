function Ad=dev(A)

I=eye(size(A)); 
Ad=A-((trace(A)./3)*I);