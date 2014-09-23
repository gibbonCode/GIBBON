function M_rice=imrician(M,s)

% function M_rice=imrician(M,s)
% ------------------------------------------------------------------------
% IMCRICIAN Random samples from the Rice/Rician probability distribution.
%
% R ~ Rice(v, s) if R = sqrt(X^2 + Y^2), where X ~ N(v*cos(a), s^2) and
% Y ~ N(v*sin(a), s^2) are independent normal distributions (any real a).
%
% Reference: http://en.wikipedia.org/wiki/Rice_distribution
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/04/2009
% ------------------------------------------------------------------------

X_GAUSS=s.*randn(size(M))+M;
Y_GAUSS=s.*randn(size(M));
M_rice = hypot(X_GAUSS,Y_GAUSS); 