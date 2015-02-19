function VI=biharm_spline_fit_ND(X,V,XI)

%   Reference:  David T. Sandwell, Biharmonic spline
%   interpolation of GEOS-3 and SEASAT altimeter
%   data, Geophysical Research Letters, 2, 139-142,
%   1987.  Describes interpolation using value or
%   gradient of value in any dimension.

%% Distance between points in x,y,z

d=dist(X,X'); 
d(1:size(d,1)+1:numel(d)) = ones(1,size(d,1)); %Replace zeros on diagonal with ones

%% Determine weights for interpolation
g = (d.^2) .* (log(d)-1);   % Green's function.
g(1:size(d,1)+1:numel(d)) = zeros(size(d,1),1); % Fixup value of Green's function along diagonal

weights = g \ V(:);

D=dist(XI,X'); L=D==0; D(L)=1; %Distance between points in x,y,z and xi,yi,zi

G=(D.^2) .* (log(D)-1);   % Green's function.
G(L)=0;

VI=G * weights;

