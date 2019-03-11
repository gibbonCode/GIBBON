function [M_median,n]=medfilt3(M,IND,k,v)

%Creating spherical mask environment with radius k
[kI,kJ,kK]=cart2im(k,k,k,v);
kI=round(kI); kJ=round(kJ); kK=round(kK);
[MASK_J,MASK_I,MASK_K]=meshgrid(-kJ:kJ,-kI:kI,-kK:kK);
% [MASK_X,MASK_Y,MASK_Z]=im2cart(MASK_I,MASK_J,MASK_K,v);
MASK_X=MASK_J.*v(2); MASK_Y=MASK_I.*v(1); MASK_Z=MASK_K.*v(3);

R=sqrt(MASK_X.^2 + MASK_Y.^2 + MASK_Z.^2);

Lv=R<=k;
MASK_I=MASK_I(Lv); MASK_J=MASK_J(Lv); MASK_K=MASK_K(Lv);

%Getting mask indices
[IND_mask]=maskfind(M,IND(:),MASK_I(:),MASK_J(:),MASK_K(:));
L_valid=IND_mask>0;
INT_valid=M(IND_mask(L_valid));
INT_mask=nan(size(IND_mask));
INT_mask(L_valid)=INT_valid;

%Calculating median, ignoring NaN's
M_median=nanmedian(INT_mask,2);

%Calculating number of elements used in median calculation
n=sum(~isnan(INT_mask),2);

%         %Creating spherical mask environment
%         k=k+iseven(k);
%         k_offset=round(k/2)-1;
%         [MASK_J,MASK_I,MASK_K]=meshgrid(-k_offset:k_offset);
%         R=sqrt(MASK_J.^2 + MASK_I.^2 + MASK_K.^2);
%         Lv=R<=k_offset;
%         MASK_I=MASK_I(Lv); MASK_J=MASK_J(Lv); MASK_K=MASK_K(Lv);
%
%         %Getting mask indices
%         [IND_mask]=maskfind(M,IND(:),MASK_I(:),MASK_J(:),MASK_K(:));
%         L_valid=IND_mask>0;
%         INT_valid=M(IND_mask(L_valid));
%         INT_mask=nan(size(IND_mask));
%         INT_mask(L_valid)=INT_valid;
%
%         %Calculating median, ignoring NaN's
%         M_median=nanmedian(INT_mask,2);
%
%         %Calculating number of elements used in median calculation
%         n=sum(~isnan(INT_mask),2);
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
