function [Vs]=smoothTaubin(F,V,nSmooth,lambda,mu)

% function [Vs]=smoothTaubin(F,V,nSmooth,lambda,mu)
% ------------------------------------------------------------------------
%
%
%
%
% ------------------------------------------------------------------------

%% Parse input

numVertices=size(V,1);

%% Get edge-vertex and vertex-vertex connectivity

E=patchEdges(F,0); %The non-unique edge set
E_sort=sort(E,2); %Sorted in column dir so 1 2 looks the same as 2 1
indEdges=sub2indn(numVertices*ones(1,2),E_sort); %Create "virtual" indices
[~,ind1,~]=unique(indEdges); %Get indices for unique edges
edgeVertexConnectivity=E(ind1,:); %Get unique edges

EV=[edgeVertexConnectivity;fliplr(edgeVertexConnectivity)]; %Non-unique edges
vertexVertexConnectivity=sparse(EV(:,1),EV(:,2),EV(:,2),numVertices,numVertices);
vertexVertexConnectivity=sort(vertexVertexConnectivity,2,'descend');
[~,J,~] = find(vertexVertexConnectivity);
vertexVertexConnectivity=full(vertexVertexConnectivity(:,1:max(J)));

logicValid=vertexVertexConnectivity>0;

%% Taubin smooth steps

Vs=V; %Allocate smoothed coordinate array
X_LAP=nan(size(vertexVertexConnectivity));
ind_X=vertexVertexConnectivity(logicValid);
ind_X_LAP=find(logicValid);
for qn=1:1:nSmooth %Loop for number of smooth iterations
    for q=1:1:size(V,2) %Loop over dimensions
        %Gaussian (Laplacian) step
        Xg=smoothCoord(V(:,q),X_LAP,ind_X_LAP,ind_X,lambda);

        %"Anti-Gaussian" step
        Xs=smoothCoord(Xg,X_LAP,ind_X_LAP,ind_X,mu);

        %Smoothed coordinates
        Vs(:,q)=Xs;
    end
    V=Vs;
end

end

function [Xs]=smoothCoord(X,X_LAP,ind_X_LAP,ind_X,p)

X_LAP(ind_X_LAP)=X(ind_X);
X_mean=mean(X_LAP,2,'omitnan');
dX=(X_mean-X);
Xs=X+p.*dX;

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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
