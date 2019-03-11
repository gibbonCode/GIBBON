function [F_fixed,L_fixed]=fixNormalsOutward(F,V,fixOpt)

%Assumes F and V describe a shape that can be appropirately defined using a
%spherical or polar coordinate system, as indicated by the choice fixOpt.
%For instance a closed spherical shape for open cylindrical shape. 
%
% The function uses the fact that the radius of the normal vector origin
% should be smaller than the radius of the normal vector tip to determine
% whether or not a normal vector needs to be flipped. 

%Derive current face normals
[N,V_starts]=patchNormal(F,V);

%Vector lengths may need to be scaled, scaling factor should be
%arbitrary if and 100 should work is shape is sufficiently
%spherical. For severly "deformed" or noisy shapes this may
%have to be smaller, but one could argue that then the
%spherical coordinate system mapping is not appropriate.
f=100;
switch fixOpt
    case 'p' %polar
        R_starts = sqrt(V_starts(:,1).^2 + V_starts(:,2).^2); %Start radii
        N=(min(R_starts(:))/f)*N; %Scaling normal lengths
        V_ends=V_starts+N; %vector end points
        R_ends = sqrt(V_ends(:,1).^2 + V_ends(:,2).^2); %End radii        
    case 's' %spherical         
        R_starts = sqrt(V_starts(:,1).^2 + V_starts(:,2).^2 + V_starts(:,3).^2); %Start radii
        N=(min(R_starts(:))/f)*N; %Scaling normal lengths
        V_ends=V_starts+N; %vector end points
        R_ends = sqrt(V_ends(:,1).^2 + V_ends(:,2).^2 + V_ends(:,3).^2); %End radii
end
L_fixed=R_ends<R_starts;
F_fixed=F;
F_fixed(L_fixed,:)=fliplr(F(L_fixed,:)); %Invert faces whose normal vector end radii are smaller than their origin

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
