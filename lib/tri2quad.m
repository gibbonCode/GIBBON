function [Fq,Vq]=tri2quad(Ft,Vt)

V1=Vt(Ft(:,1),:);
V2=Vt(Ft(:,2),:);
V3=Vt(Ft(:,3),:);
V4=(V1+V2)./2;
V5=(V2+V3)./2;
V6=(V3+V1)./2;
V7=(V1+V2+V3)./3;

E=[Ft reshape([size(Vt,1)+1:size(Vt,1)+4*size(Ft,1)]',size(Ft,1),4)];

Vq=[Vt;V4;V5;V6;V7];
Fq=[E(:,1) E(:,4) E(:,7) E(:,6);...
    E(:,4) E(:,2) E(:,5) E(:,7);...
    E(:,6) E(:,7) E(:,5) E(:,3)];

%Removing double vertices
[~,ind1,ind2]=unique(pround(Vq,5),'rows');
Vq=Vq(ind1,:);
Fq=ind2(Fq);




