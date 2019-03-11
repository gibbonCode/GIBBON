function [I,J,K]=MRcart2imNeg(Xc,Yc,Zc,v,OR,rn,cn)

%--------------------------------------------------------------------------
% function [I,J,K]=MRcart2im(X,Y,Z,v,OR,r,c)
% 
% This function calculates the image coordinates I,J,K of the cartesian 
% coordinates X,Y,Z based on the voxel dimensions v, the origin OR (cartesian
% coorindates of the point I=1, J=1, K=1), the offsets r anc c.
%
% v=N_info.PixelSpacing;
% v(3)=N_info.SliceThickness;
% OR=N_info.ImagePositionPatient;
% r=N_info.ImageOrientationPatient(4:6);
% c=N_info.ImageOrientationPatient(1:3);
%
% Based on:
% http://cmic.cs.ucl.ac.uk/fileadmin/cmic/Documents/DavidAtkinson/DICOM_6up.pdf
% and the file "TranTransform matrix between two dicom image coordinates"
% from the matlab central file exchange written by Alper Yaman. 
%-------------------------------------------------------------------------- 

reForm=~isvector(Xc);

%Switch convention
X=Yc; Y=Xc; Z=Zc;

size_M=size(Xc);

% r => row direction vector
r=vecnormalize(cn);

% c => column direction vector
c=vecnormalize(rn);

%Determine s => slice direction vector N.B. defined as cross(c,r) not
%cross(r,c) 
s=cross(c',r'); %Determine s => slice direction vector
s=vecnormalize(s);

X=X(:)';Y=Y(:)';Z=Z(:)';

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

XYZ=ones(4,length(X));
XYZ(1,:)=X;
XYZ(2,:)=Y;
XYZ(3,:)=Z;

IJK=(M\XYZ)';

I=IJK(:,1)+1; J=IJK(:,2)+1; K=IJK(:,3)+1;

if reForm %If the input is not a vector reshape
    I=reshape(I,size_M);
    J=reshape(J,size_M);
    K=reshape(K,size_M);
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
