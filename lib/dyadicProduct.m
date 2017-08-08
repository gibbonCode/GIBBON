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
% ********** _license boilerplate_ **********
% 
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
