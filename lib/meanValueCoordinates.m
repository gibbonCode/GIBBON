function [C]=meanValueCoordinates(V_cage,V_interp)

% function [C]=meanValueCoordinates(V_cage,V_interp)
% ------------------------------------------------------------------------
% WORK IN PROGRESS
% 
% Based on: https://www.numerical-tours.com/matlab/meshdeform_4_barycentric/
% 
% 2023/04/26
% ------------------------------------------------------------------------

%% Parse input 

%% Derived parameters

toleranceLevel=eps(0);
nInterp=size(V_interp,1); 
nCage=size(V_cage,1); 

%% Process mean value coordinates
C = zeros(nInterp,nCage);
D = nan(nInterp,nCage);
for qCagePoint=1:nCage
    V_now = V_cage(qCagePoint,:);
    U_mesh = V_now(ones(nInterp,1),:)-V_interp;
    nb = vecnormalize(U_mesh);
    nb(:,3)=0;    
    d = sqrt(sum(U_mesh.^2,2));    
    D(:,qCagePoint)=d;    
    
    s = 1; %Initialise as 1
    for j=mod([qCagePoint-2,qCagePoint],nCage)+1
        vj = V_cage(j,:);
        na = vecnormalize( repmat(vj,[nInterp 1])-V_interp );
        na(:,3)=0;
        
        % sign
        si = s*sign(cross(na,nb,2));
        si=si(:,3);
        
        % angle
        dp = dot(na,nb,2);
        theta = si .* acos(rescale(dp,-1,1));
        
        % add tangent of half angle
        C(:,qCagePoint) = C(:,qCagePoint) + reshape( tan(theta/2) ./ d, [nInterp 1]);

        s = -s; %Flip sign
    end
end

logicReplace=D<toleranceLevel;
C(any(logicReplace,2),:)=0;
C(logicReplace)=1;

C = C ./ sum(C,2); % Normalization

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
