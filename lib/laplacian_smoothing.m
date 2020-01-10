function V=laplacian_smoothing(V,IND_V,L,n)

[I,J,v] = find(IND_V);

for q=1:n
    Xp=accumarray({I,J},V(v,1),size(IND_V),[],NaN);
    Yp=accumarray({I,J},V(v,2),size(IND_V),[],NaN);
    Zp=accumarray({I,J},V(v,3),size(IND_V),[],NaN);
    Vp=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)];
    V=V+L.*(Vp-V);
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
