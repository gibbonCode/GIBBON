function [Fm,Vm]=cap_patchcylinder(F1,V1,F2,V2,nr,nz)

% function [Fm,Vm]=cap_patchcylinder(F1,V1,F2,V2,nr,nz)
% ------------------------------------------------------------------------
% This function assumes the inputs F1, V1 and F2, V2 define the faces and
% vertices of two cylinders and closes to top and bottom faces by
% connecting the cylinders together. The input nr defines the number
% of radial points to add for the connection. The input nz defines the
% number of steps used in the z-direction for the cylinders. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

Fm=[F1; (F2+size(V1,1))]; 
Vm=[V1; V2];

%Coordinates created from nr*nz matrix 

%Top edge indices
It=ones(nr-1,1); Jt=(1:1:nr-1)';
INDt=sub2ind([nz,nr-1],It,Jt);
% Vt1=V1(INDt,:); Vt2=V2(INDt,:); Vt=[Vt1; Vt2];

%Bottom edge indices
Ib=nz.*ones(nr-1,1); Jb=(1:1:nr-1)';
INDb=sub2ind([nz,nr-1],Ib,Jb);
% Vb1=V1(INDb,:); Vb2=V2(INDb,:); Vb=[Vb1; Vb2];

Fp=(ones(nr-1,1)*[1 2 2 1])+(((1:1:nr-1*ones)'-1)*ones(1,4)); Fp(Fp==nr)=1; %Quad order

Fq=[INDt(Fp(:,1:2)) INDt(Fp(:,3:4))+size(V1,1)];
Ftt=[Fq(:,[1 3 2]); Fq(:,[1 4 3])]; %Tri order

Fq=[INDb(Fp(:,1:2)) INDb(Fp(:,3:4))+size(V1,1)];
Ftb=[Fq(:,[1 3 2]); Fq(:,[1 4 3])]; %Tri order

Fm=[Fm; Ftt; Ftb];

end