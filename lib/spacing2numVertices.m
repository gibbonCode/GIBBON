function nP=spacing2numVertices(F,V,pointSpacing)

A=sum(patchArea(F,V)); %Total area
At=(pointSpacing.^2*sqrt(3))/4; %Theoretical area of equilateral triangle
NF=(A/At);

E=patchEdges(F,1); %Edges
nE=size(E,1); %Number of edges
nF=size(F,1); %Number of faces
nV=size(V,1); %Number of vertices

%Compute Euler characteristic 
X=size(V,1)-size(E,1)+size(F,1); 

nRefScalar=(log(NF)-log(nF))/log(4);

if nRefScalar<0
    nRef=floor(nRefScalar);
    nRange=0:-1:nRef;
else
    nRef=ceil(nRefScalar);
    nRange=0:1:nRef;
end

nvR=nV.*ones(numel(nRange),1);
neR=nE.*ones(numel(nRange),1);

nfR=nF*4.^nRange';
for q=2:1:numel(nRange)
    if nRef>0                
        nvR(q)=nvR(q-1)+neR(q-1);        
    elseif nRef<0
        nvR(q)=(nvR(q-1)+X-nfR(q))/2;
    end    
    neR(q)=-X+nfR(q)+nvR(q);
end

% [nvR -neR nfR]
% sum([nvR -neR nfR],2)

nP=ceil(interp1(nfR(:),nvR(:),NF,'pchip'));
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
