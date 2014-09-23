function P=ricianpdf(A,M,S)

P= (((M./(S.^2)).*exp((-1.*((A.^2)+(M.^2)))./(2.*(S^2)))).*(besseli(0, ((A.*M)./(S.^2))))).*(M>=0);