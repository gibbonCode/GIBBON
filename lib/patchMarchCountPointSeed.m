function [indSeed,marchCount,seedIndex]=patchMarchCountPointSeed(F,V,indStart,n)

[~,IND_V]=tesIND(F,V,0);
optStruct.IND_V=IND_V;

numStarts=numel(indStart);
indSeed=nan(n,1);
indSeed(1:numStarts)=indStart; 

for q=numStarts:1:n+1
    [marchCount,seedIndex]=patchMarchCount(F,V,indSeed(~isnan(indSeed)),optStruct);
    [~,indAdd]=max(marchCount);
    if q<=n
        indSeed(q)=indAdd;
    end
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
