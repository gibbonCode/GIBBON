function C=IJKL_JK_doubleContraction(A,B)

% clear all; close all; clc;
% A=rand(3,3,3,3);
% B=rand(3,3);

sizA=size(A);
sizB=size(B);

siz4=3.*ones(1,4);
siz2=3.*ones(1,2);

if ndims(A)==4 && ismatrix(B) %4d array and 2d array
    if all(sizA==siz4) && all(sizB==siz2) %4th order tensor and 2nd order tensor
        
        %Creating indices for tensors
        IJKL=gcombvec(1:3,1:3,1:3,1:3)';
        I=IJKL(:,1); J=IJKL(:,2); K=IJKL(:,3); L=IJKL(:,4);
        IND_IJKL=sub2ind(siz4,I,J,K,L); %Linear indices for 4th order tensor A
        IND_JK=sub2ind(siz2,J,K); %Linear indices for 2nd order tensor B
        
        if isa(A,'double') && isa(B,'double') %if both are doubles
            C=zeros(siz2);
        elseif isa(A,'sym') || isa(B,'sym') %if one of the is a symbolic
            C=sym(zeros(siz2));
        else
            %TO DO! Add warning here
        end
        
        %Create basis vectors and dyadic product N.B. these are chosen here!
        E=eye(3,3); 
        
        %Sum across indices IJKL
        for i=1:1:numel(IND_IJKL)
            EIL=(E(I(i),:).')*E(L(i),:); %dyadic product of basis vectors eI and eL
            c=A(IND_IJKL(i)).*B(IND_JK(i)).*EIL; %sub-c
            C=C+c; %summed C
        end
    end
else
    %TO DO! Add warning here
end


 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
