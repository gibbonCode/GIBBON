function VI=biharmonicSplineInterpolation(X,V,XI)

% function VI=biharmonicSplineInterpolation(X,V,XI)
% ------------------------------------------------------------------------
%
% The |biharmonicSplineInterpolation| function is an expansion to
% n-dimensions and scattered data of the biharmonic spline interpolation
% method of the |griddata| function (method 'v4').
%
%   Reference:  David T. Sandwell, Biharmonic spline interpolation of
%   GEOS-3 and SEASAT altimeter data, Geophysical Research Letters, 2,
%   139-142, 1987.  
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

%%
% Convert row to column data

if isrow(XI)
   XI=XI'; 
end

if isrow(X)
    X=X'; 
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

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
