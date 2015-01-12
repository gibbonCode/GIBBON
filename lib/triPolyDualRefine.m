function [Ft,Vt,Ct,indIni]=triPolyDualRefine(F,V)

[Vd,~,Fds]=patch_dual(V,F);

Cd=(1:1:size(Fds,1))';
Cds=Cd(:,ones(size(Fds,2),1));

Fds_t=Fds'; 
Cds_t=Cds';

Vt=[Vd;V];

Q1=Fds_t(:);
L1=Q1>0; 
Q1=Q1(L1);
Ct=Cds_t(:);
Ct=Ct(L1);

Q2=Fds(:,[2:size(Fds,2) 1])';
Q2=Q2(:);
L2=Q2>0; 
Q2=Q2(L2);

q=1:size(Fds,1);
Q3=q(ones(size(Fds,2),1),:)+size(Vd,1);
Q3=Q3(:);
A=Q3;
Q3=Q3(L1);
Ft=[Q1 Q2 Q3];

indIni=(size(Vd,1)+1):(size(Vd,1)+size(V,1));

