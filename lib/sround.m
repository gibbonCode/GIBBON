function Xr=sround(X,NumSigDigits)

logicSign=X<0;
X=abs(X);
scaleFactors=10.^(floor(log10(X))+1);
Xr=pround(X./scaleFactors,NumSigDigits).*scaleFactors;
Xr(logicSign)=-Xr(logicSign); clc

end