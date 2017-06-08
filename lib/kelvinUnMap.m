function C=kelvinUnMap(cKelvin)

x=[1 1 1 sqrt(2) sqrt(2) sqrt(2)]'; %conversion Voigt to Kelvin form

if isvector(cKelvin) 
   cVoigt=cKelvin(:)./x;
else
    X=x.'*x;
    cVoigt=cKelvin./X;
end

C=voigtUnMap(cVoigt);
