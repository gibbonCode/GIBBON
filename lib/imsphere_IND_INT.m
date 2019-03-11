function [IND INT]=imsphere_IND_INT(r,v_low,res,size_IM,Xp,Yp,Zp)

% function [I J K]=imsphere_IND_INT(r,v_low,res,size_IM,Xp,Yp,Zp)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 07/08/2008
% ------------------------------------------------------------------------

%% Setting up voxel dimensions and voxel radius dimensions
v_high=v_low/res;
r_v_high=r/v_high(1);
r_v_low=r/v_low(1);

%% Calculating image coordinate of point in high resolution image
[Ip_high Jp_high Kp_high]=cart2im(Xp,Yp,Zp,v_high);

%% Finding the index of the voxel in which point is found or is closest to
Ir_high=round(Ip_high); Jr_high=round(Jp_high); Kr_high=round(Kp_high); 

%% Calculating shift from the centre of this voxel
Ic_high=Ir_high-Ip_high; Jc_high=Jr_high-Jp_high; Kc_high=Kr_high-Kp_high;

%% Finding indices of points in a sphere at specified location
[I_sph_high J_sph_high K_sph_high]=sphere_index(r_v_high,Ic_high,Jc_high,Kc_high);

%Calculating image coordinate of point in low resolution image
[Ip_low Jp_low Kp_low]=cart2im(Xp,Yp,Zp,v_low);
%Finding the index of the voxel in which point is found
Ir_low=round(Ip_low); Jr_low=round(Jp_low); Kr_low=round(Kp_low); 
%Calculating shift from the centre of this voxel
Ic_low=Ir_low-Ip_low; Jc_low=Jr_low-Jp_low; Kc_low=Kr_low-Kp_low;

%Calculating shift
I_low_high_shift=((Ir_low*res)-(res/2)+0.5)-Ir_high;
J_low_high_shift=((Jr_low*res)-(res/2)+0.5)-Jr_high;
K_low_high_shift=((Kr_low*res)-(res/2)+0.5)-Kr_high;

size_M=(res+(2*round(r_v_low+1)*res))*ones(1,3);
M=ones(size_M);
IJK_middle=round(size_M/2);

IJK_voxel_center(1)=IJK_middle(1)-I_low_high_shift;
IJK_voxel_center(2)=IJK_middle(2)-J_low_high_shift;
IJK_voxel_center(3)=IJK_middle(3)-K_low_high_shift;
%round_IJK_voxel_center=round(IJK_voxel_center);

I_sph_high=round(I_sph_high+IJK_voxel_center(1)); 
J_sph_high=round(J_sph_high+IJK_voxel_center(2)); 
K_sph_high=round(K_sph_high+IJK_voxel_center(3));
IND = sub2ind(size_M,I_sph_high,J_sph_high,K_sph_high);
M(IND)=0;

%disp(['High resolution sphere simulated using ', num2str(numel(IND)), ' voxels']);

%Calculating low resolution image
[M]=voxelate(M,res);

IND=find(M<1); 
INT=M(IND);
[I,J,K]=ind2sub(size(M),IND); 
IJK_middle=round(size(M)/2);
I=I-IJK_middle(1); J=J-IJK_middle(2); K=K-IJK_middle(3);

I=I+Ir_low;
J=J+Jr_low;
K=K+Kr_low;

[IND]=sub2ind(size_IM,I,J,K); 

%% END
 
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
