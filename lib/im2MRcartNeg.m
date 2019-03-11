function [Xc,Yc,Zc]=im2MRcartNeg(I,J,K,v,OR,rn,cn)

%--------------------------------------------------------------------------
% function [X,Y,Z]=im2MRcart(I,J,K,v,OR,r,c)
% 
% This function calculates the cartesian coordinates X,Y,Z of the image
% coordinates I,J,K based on the voxel dimensions v, the origin OR (cartesian
% coorindates of the point I=1, J=1, K=1), the offsets r anc c.
%
% v=N_info.PixelSpacing;
% v(3)=N_info.SliceThickness;

% ORIGIN:
% The x, y, and z coordinates of the upper left hand corner (center of the
% first voxel transmitted) of the image, in mm. 
% OR=N_info.ImagePositionPatient; 

% ROW DIRECTION
% r=N_info.ImageOrientationPatient(4:6);

% COLUMN DIRECTION
% c=N_info.ImageOrientationPatient(1:3);
%
% Based on:
% http://cmic.cs.ucl.ac.uk/fileadmin/cmic/Documents/DavidAtkinson/DICOM_6up.pdf
% and the file "TranTransform matrix between two dicom image coordinates"
% from the matlab central file exchange written by Alper Yaman. 
%-------------------------------------------------------------------------- 

reForm=~isvector(I); 

size_M=size(I);

% r => row direction vector
r=vecnormalize(cn);

% c => column direction vector
c=vecnormalize(rn);

%Determine s => slice direction 
s=cross(c',r'); %Determine s => slice direction vector
s=vecnormalize(s);

I=(I(:)-1)';J=(J(:)-1)';K=(K(:)-1)'; %Shift so first voxel is at 0 0 0

%Translation
T  = [1 0 0 OR(1);...
      0 1 0 OR(2);...
      0 0 1 OR(3);...
      0 0 0 1]; 
%Rotation
R  = [r(1) c(1) s(1) 0;...
      r(2) c(2) s(2) 0;...
      r(3) c(3) s(3) 0;...
      0    0    0    1];
%Scaling  
S  = [v(1) 0    0    0;...
      0    v(2) 0    0;...
      0    0    v(3) 0;...
      0    0    0    1];
  
T0 = eye(size(T));

M  = T * R * S * T0; %The transformation matrix

IJK=ones(4,length(I));
IJK(1,:)=I;
IJK(2,:)=J;
IJK(3,:)=K;
IJK=double(IJK);

XYZ=(M*IJK)';

X=XYZ(:,1); Y=XYZ(:,2); Z=XYZ(:,3);

if reForm %If the input is not a vector reshape
    X=reshape(X,size_M);
    Y=reshape(Y,size_M);
    Z=reshape(Z,size_M);
end

%Switch convention
Xc=Y; Yc=X; Zc=Z;


 
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
