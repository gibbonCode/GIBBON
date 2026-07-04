function [D1]=triSurfSetDist(F1,V1,F2,V2,distMetric)

% function [D1]=triSurfSetDist(F1,V1,F2,V2,distMetric)
% ------------------------------------------------------------------------
%
% Change log: 
% 2021/07/22 KMM Removed waitbar for ray-tracing method
% ------------------------------------------------------------------------

%%

switch distMetric
    case 'dist'
        D1=minDist(V1,V2);
    case 'ray'
        [~,~,N1]=patchNormal(F1,V1);
        
        %Ray trace 1 onto 2        
        optionStructRayTrace.tolEps = 1e-6;
        optionStructRayTrace.triSide = 0;
        optionStructRayTrace.rayType = 'ray'; % or line
        optionStructRayTrace.exclusionType = 'normal'; % or exclusive
        optionStructRayTrace.paired=0; % 1 for paired, 0 for non-paired

        [~,indIntersect,d1_all]=triSurfRayTrace(V1,N1,F2,V2,optionStructRayTrace);
  
        numIntersect=size(indIntersect,1);

        S=sparse(indIntersect(:,1),(1:numIntersect)',1,size(V1,1),numIntersect,numIntersect);
        d1_sparse=sparse(indIntersect(:,1),(1:numIntersect)',abs(d1_all),size(V1,1),numIntersect,numIntersect);
        
        D1=full(spmin(d1_sparse,[],2,'includenan',S~=0,1));        
        
    case 'dist-ray'
        [D1d]=triSurfSetDist(F1,V1,F2,V2,'dist');
        [D1r]=triSurfSetDist(F1,V1,F2,V2,'ray');
        D1=gnanmin([D1d D1r],[],2);     
    case 'near-norm'
        [~,indMin]=minDist(V1,V2);
        [~,~,Nv]=patchNormal(F2,V2);
        N=Nv(indMin,:);
        W=V1-V2(indMin,:);
        D1=abs(dot(N,W,2));
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
