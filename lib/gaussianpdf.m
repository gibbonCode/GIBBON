function P=gaussianpdf(A,M,S)

P= ((1./(S.*sqrt(2*pi))).*exp(-1.*((M-sqrt((A^2)+(S^2))).^2./(2.*(S.^2))))).*(M>=0);