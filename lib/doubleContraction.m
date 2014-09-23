function C=doubleContraction(A,B)
% clear all; close all; clc;
% A=rand(3,3,3,3);
% B=rand(3,3);

sizA=size(A);
sizB=size(B);

siz4=3.*ones(1,4);
siz2=3.*ones(1,2);

if ismatrix(A) && ismatrix(B) %Two 2D arrays
    if all(sizA==siz2) && all(sizB==siz2) %Two 2nd order tensors
        C=trace(A*(B.')); %C=Aij.*Bij
    end
elseif ndims(A)==4 && ismatrix(B) %4d array and 2d array
    if all(sizA==siz4) && all(sizB==siz2) %4th order tensor and 2nd order tensor
        
        %Creating indices for tensors
        IJKL=combvec(1:3,1:3,1:3,1:3)';
        I=IJKL(:,1); J=IJKL(:,2); K=IJKL(:,3); L=IJKL(:,4);
        IND_IJKL=sub2ind(siz4,I,J,K,L); %Linear indices for 4th order tensor A
        IND_KL=sub2ind(siz2,K,L); %Linear indices for 2nd order tensor B
        
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
            EIJ=(E(I(i),:).')*E(J(i),:); %dyadic product of basis vectors eI and eJ
            c=A(IND_IJKL(i)).*B(IND_KL(i)).*EIJ; %sub-c
            C=C+c; %summed C
        end
    end
else
    %TO DO! Add warning here
end

end

