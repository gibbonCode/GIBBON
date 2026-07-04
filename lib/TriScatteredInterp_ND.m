function [DI]=TriScatteredInterp_ND(P,U,Pi,interpMethod)

% function [DI]=TriScatteredInterp_ND(P,U,Pi,interpMethod)
% ------------------------------------------------------------------------
%
%
%
% ------------------------------------------------------------------------

%%

%Deal with older Delaunay type input
if isa(P,'DelaunayTri') || isa(P,'delaunayTriangulation') %if DT is not a delaunay tesselation
    P=P.Points; %Override as just the points
end

%Process interpolation for each dimension
DI=nan(size(Pi,1),size(U,2)); %Allocate DI
for q=1:size(U,2)% loop over dimensions
    switch interpMethod
        case 'nat_near' %natural in chull, neirest outside chull
            [DI(:,q),~]=TriScatteredInterp_nat_near(P,U(:,q),Pi);
        otherwise %TriScatterInterp can handle other methods
            F = scatteredInterpolant(P,U(:,q),interpMethod); %Construct interpolator
            DI(:,q)=F(Pi); %Interpolate
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
