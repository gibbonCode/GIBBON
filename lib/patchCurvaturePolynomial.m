function [U1,U2,K1,K2,H,G] = patchCurvaturePolynomial(F,V)

% function [U1,U2,K1,K2,H,G] = patchCurvaturePolynomial(F,V)
% ------------------------------------------------------------------------
%
% This function computes the mesh curvature at each vertex for the input mesh 
% defined by the face `F` and the vertices `V`. A local polynomial is fitted to 
% each point's "Laplacian umbrella" (point neighbourhood), and the curvature of 
% this fitted form is derived. 
% 
% The reference below [1] provides more detail on the algorithm. In addition, this 
% implementation was created with the help of this helpful document:
% https://github.com/alecjacobson/geometry-processing-curvature/blob/master/README.md, 
% which features a nice overview of the theory/steps involved in this algorithm. 
% 
% References 
% 1. [F. Cazals and M. Pouget, _Estimating differential quantities using polynomial fitting of osculating jets_, Computer Aided Geometric Design, vol. 22, no. 2, pp. 121-146, Feb. 2005, doi: 10.1016/j.cagd.2004.09.004](https://doi.org/10.1016/j.cagd.2004.09.004)
% 
% ------------------------------------------------------------------------

%%
m = size(V,1);
CC = patchConnectivity(F,V,{'vv'});
con_V2V = CC.vertex.vertex;
[~,~,NV]=patchNormal(F,V);
nz = [0.0 0.0 1.0]; % z-vector

K1 = zeros(m,1);
K2 = zeros(m,1);
U1 = zeros(m,3);
U2 = zeros(m,3);
for q = 1:1:m
    n = NV(q,:);
    [a,d]=vectorOrthogonalPair(n);
    Q=[a; d; n];
    ind = con_V2V(q,:); ind = ind(ind>0);
    vr = (V(ind,:)-V(q,:))*Q';
    
    % Set up polynomial fit
    T = zeros(length(ind),5);
    w = zeros(length(ind),1);
    for i = 1:length(ind)
        T(i,:) = [vr(i,1),vr(i,2),vr(i,1)^2,vr(i,1)*vr(i,2),vr(i,2)^2];
        w(i) = vr(i,3);
    end
    a = T\w;

    E = 1.0 + a(1)^2;
    F = a(1)*a(2);
    G = 1.0 + a(2)^2;
    d = sqrt(a(1)^2+1.0+a(2)^2);
    e = (2.0*a(3)) / d;
    f =       a(4) / d;
    g = (2.0*a(5)) / d;

    S = -[e f; f g] * inv([E F; F G]);
    [u,k] = eig(S); % Eigen decomposition to get first/second eigenvalue and vectors    

    % Store derived quantities
    K1(q) = k(2,2);
    K2(q) = k(1,1);
    U1(q,:) = [u(1,2) u(2,2) 0.0]*Q;
    U2(q,:) = [u(1,1) u(2,1) 0.0]*Q;   

end

H = 0.5 * (K1+K2); % Mean curvature
G = K1.*K2; % Gaussian curvature

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
