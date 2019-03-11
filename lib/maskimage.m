function [mask_result]=maskimage(M,IND,I,J,K,W)

% function [M_mask]=maskimage(M,IND,I,J,K,W);
% ------------------------------------------------------------------------
% This function applies the mask defined by the centered indices I,J,K and
% the weights W to the image M and returns the masked image in the vector
% mask_result
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/09/2008
% ------------------------------------------------------------------------

%% Padding array

%Extra space is added such that the boundary elements can also be masked
repSize=round([max(I)-min(I) max(J)-min(J) max(K)-min(K)]/2)+1;
[M_pad]=padrep(M,repSize);

%% Setting up masking matrices

%Converting the index numbers IND of M to the corresponding index number
%IND_pad of M_pad
size_M_pad=size(M_pad);
[I_pad J_pad K_pad]=ind2sub(size(M),IND);
I_pad=I_pad+repSize(1);
J_pad=J_pad+repSize(2);
K_pad=K_pad+repSize(3);
IND_pad=sub2ind(size_M_pad,I_pad,J_pad,K_pad);

%Converting the centered indices I,J,K to centered linear indEx numbers
%L_ind
IJK_M=round(size_M_pad/2);
I_M=I+IJK_M(1); J_M=J+IJK_M(2); K_M=K+IJK_M(3);
n_M = sub2ind(size_M_pad,I_M,J_M,K_M);
L_ind = n_M-sub2ind(size_M_pad,IJK_M(1),IJK_M(2),IJK_M(3));

%Setting up masking matrix, which initially numel(IND) x numel(W) in
%size, each row corresponds to a voxel number and each column corresponds
%to the number of the voxels used to define the masked value for this voxel
mask_voxels=ones(length(IND_pad),1)*L_ind'+IND_pad*ones(1,size(L_ind,1));

%Defining weights matrix. This matrix has the same size as the masking
%matrix but contains the weights for each voxel on the columns
if size(W,1)==1
    mask_weights=(ones(length(IND_pad),1)*W);
else
    mask_weights=(W*ones(1,length(IND_pad)))';
end

%% Appling mask to image

mask_data=M_pad(mask_voxels).*mask_weights;
mask_result=sum(mask_data,2);

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
