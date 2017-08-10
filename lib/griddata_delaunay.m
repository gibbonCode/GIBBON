function [zi] = griddata_delaunay(x,y,z,xi,yi,DT,method)

%Simple hack of griddata. Error checking removed, delaunay triangulation
%provided. 

%%

tri=DT.Triangulation;

%% Checking method
switch lower(method),
    case 'linear'
        zi = linear(x,y,z,xi,yi,tri);
    case 'cubic'
        if(isreal(z))
            zi = cubic(x,y,z,xi,yi,tri);
        else
            zi = complex(cubic(x,y,real(z),xi,yi,tri),cubic(x,y,imag(z),xi,yi,tri));
        end
    case 'nearest'
        zi = nearest(x,y,z,xi,yi,tri);
    otherwise
        error('MATLAB:griddata:UnknownMethod', 'Unknown method.');
end


%% LINEAR
function zi = linear(x,y,z,xi,yi,tri)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

% Only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);

% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
    (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
    (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
    (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
    (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);

z = z(:).'; % Treat z as a row so that code below involving
% z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

if ~isempty(out), zi(out) = NaN; end
%------------------------------------------------------------

%% CUBIC

function zi = cubic(x,y,z,xi,yi,tri)
%TRIANGLE Triangle-based cubic interpolation

%   Reference: T. Y. Yang, "Finite Element Structural Analysis",
%   Prentice Hall, 1986.  pp. 446-449.
%
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);
zi = cubicmx(x,y,z,xi,yi,tri,t);
%------------------------------------------------------------

%% NEAREST NEIGHBOR

function zi = nearest(x,y,z,xi,yi,tri)
%NEAREST Triangle-based nearest neighbor interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these a columns
x = x(:); y = y(:); z = z(:); % Treat these as columns

% Find the nearest vertex
k = dsearch(x,y,tri,xi,yi);
% k = nearestNeighbor(DT, QX)

zi = k;
d = find(isfinite(k));
zi(d) = z(k(d));
zi = reshape(zi,siz);
%----------------------------------------------------------

 
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
