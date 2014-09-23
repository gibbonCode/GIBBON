function P=rayleighpdf(M,S)

P= (((M./(S.^2)).*exp((-1.*(M.^2))./(2.*S^2)))).*(M>=0);