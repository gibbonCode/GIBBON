function [e]=erfci(a)
% The complex extension to the complimentary error function
% MATLAB does not support erfc(a) if a is a complex number. 
% erfc(a) = 1-erf(a) for a real number a. 
% erfci(a) = 1+1i*erfi(1i*a) for a complex number a. 

e = 1+1i*erfi(1i*a); 

end