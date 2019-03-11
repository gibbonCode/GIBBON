function [Vs]=triSurfSmoothFourier(F,V,numFreq)

%Based on: 
% https://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/meshproc_4_fourier/


%%
numPoints = size(V,1); %Number of points

% The combinatorial laplacian is a linear operator (thus a NxN matrix where N is the number of vertices). It depends only on the connectivity of the mesh, thus on face only.
% Compute edge list.
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])]; %Edges array
numEdges = size(E,1); %Number of edges

% Compute the adjacency matrix.
W = sparse(E(:,1),E(:,2),ones(numEdges,1),numPoints,numPoints,numEdges);
W = max(W,W');

% Compute the combinatorial Laplacian, stored as a sparse matrix.
D = spdiags(sum(W,1)',0,numPoints,numPoints);
L = D-W;

opts.disp = 0;
[U,S] = eigs(L,numFreq,'SM',opts);
S = diag(S);

% Order the eigenvector by increasing frequencies.
[~,indSort] = sort(S,'ascend');
U = real(U(:,indSort));

% Linear Approximation over the Fourier Domain is obtained by keeping only
% the low frequency coefficient. This corresponds to a low pass filtering,
% since high frequency coefficients are removed.

% Compute the projection of each coordinate vertex(i,:) on the small set of nb frequencies.
vertexSpectrum = V'*U;

% Reconstruct the mesh.
Vs = (vertexSpectrum*U')';

end
 
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
