function VI=biharmonicSplineInterpolation(X,V,XI)

% function VI=biharmonicSplineInterpolation(X,V,XI)
% ------------------------------------------------------------------------
%
% The |biharmonicSplineInterpolation| function is an expansion to
% n-dimensions and scattered data of the biharmonic spline interpolation
% method of the |griddata| function (method 'v4').
% 
%
%   Reference:  David T. Sandwell, Biharmonic spline interpolation of
%   GEOS-3 and SEASAT altimeter data, Geophysical Research Letters, 2,
%   139-142, 1987.  Describes interpolation using value or gradient of
%   value in any dimension.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

%%

if isvector(XI) && size(XI,1)==1
   XI=XI(:); 
end

if isvector(X) && size(X,1)==1
   X=X(:); 
end

%% Distances from all points in X to all points in X

D=dist(X,X'); 
D(1:size(D,1)+1:numel(D)) = ones(1,size(D,1)); % Replace zeros on diagonal with ones

%% Determine weights for interpolation
g = (D.^2) .* (log(D)-1);   % Green's function.
g(1:size(D,1)+1:numel(D)) = zeros(size(D,1),1); % Fixup value of Green's function along diagonal

W = g \ V(:);

D=dist(XI,X'); % Distance between points in X and XI
L=(D==0); D(L)=1; % Replace zeros with ones

G=(D.^2).*(log(D)-1); % Green's function.
G(L)=0;

VI=G * W;

