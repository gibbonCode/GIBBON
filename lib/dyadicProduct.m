function C=dyadicProduct(A,B,r)

%Computes dyadic products of two second order tensors yielding a
%fourth-order tensor. See Curnier et al. 1994 and Ateshian et al. 2009

%Setting up the Fourth order tensor
siz4=3.*ones(1,4);
siz2=3.*ones(1,2);

if isa(A+B,'double')
    C=zeros(siz4);
elseif isa(A+B,'sym')
    C=sym(zeros(siz4));
else
    %TO DO! Add warning here
end

%Creating indices for tensors
IJKL=gcombvec(1:3,1:3,1:3,1:3)';
I=IJKL(:,1); J=IJKL(:,2); K=IJKL(:,3); L=IJKL(:,4);
IND_IJKL=sub2ind(siz4,I,J,K,L);

IND_II=sub2ind(siz2,I,I);
IND_JJ=sub2ind(siz2,J,J);

IND_IJ=sub2ind(siz2,I,J);
IND_KL=sub2ind(siz2,K,L);

IND_IK=sub2ind(siz2,I,K);
IND_JL=sub2ind(siz2,J,L);
IND_LJ=sub2ind(siz2,L,J);

IND_IL=sub2ind(siz2,I,L);
IND_JK=sub2ind(siz2,J,K);
IND_KJ=sub2ind(siz2,K,J);

switch r
    case 1
        C(IND_IJKL)=A(IND_IJ).*B(IND_KL);
    case 2
        C(IND_IJKL)=A(IND_IK).*B(IND_JL);
    case 3
        C(IND_IJKL)=0.5.*( (A(IND_IK).*B(IND_JL)) + (A(IND_IL).*B(IND_JK)) );
    case 4
        C(IND_IJKL)=0.5.*( (A(IND_IK).*B(IND_LJ)) + (A(IND_IL).*B(IND_KJ)) );
    case 5
        C(IND_IJKL)=A(IND_IJ).*B(IND_IJ);
    case 6
        C(IND_IJKL)=A(IND_II).*B(IND_JJ);
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
