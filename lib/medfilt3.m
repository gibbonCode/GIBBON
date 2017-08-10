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



 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
