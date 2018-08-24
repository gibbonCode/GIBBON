function [VV]=triSurfVolume(F,V)

% function [VV]=triSurfVolume(F,V)
% ------------------------------------------------------------------------
%Volume derivation based on Gauss divergence theorem
%
%%% EXAMPLE
%
% r=3; %sphere radius
% n=2; %Refinements   
% [F,V,~]=geoSphere(4,r);
% 
% VVt=(4/3)*pi*r.^3; %Theoretical volume of the sphere
% [VV]=triSurfVolume(F,V); %estimate based on triangulated surface
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 26/11/2013
% 2018/08/24 Fixed bug in relation to surface not centred on centroid and
% surfaces with inhomogeneous node distributions
%------------------------------------------------------------------------
%%

% Centre on centroid
Vm=triSurfCentroid(F,V);
V=V-Vm(ones(size(V,1),1),:);

% meanV=mean(V,1);
% V=V-meanV(ones(size(V,1),1),:);

[N,~]=trinorm(F,V); %Face normals
aa=tri_area(F,V); %Areas
[Vm]=patchCentre(F,V);
Zm=Vm(:,3); %Mean Z
Nz = N(:,3); %Z component of normal
vv = aa.*Zm.*Nz; %Contributions
VV = sum(abs(vv)); %Total volume
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
