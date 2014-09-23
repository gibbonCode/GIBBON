function cKelvin=kelvinMap(c)

cVoigt=voigtMap(c);
x=[1 1 1 sqrt(2) sqrt(2) sqrt(2)]; %conversion Voigt to Kelvin form

if ~isvector(cVoigt) %assume that c is a 4th order tensor
    X=x.'*x;
    cKelvin=X.*cVoigt;
elseif isvector(cVoigt) %assume that c is a 2nd order tensor
    cKelvin=cVoigt.*x';
end
