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

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
