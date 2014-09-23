function [Lk]=logicErodeDilate(L,k,erodeDilateOption)

%Normalise structural kernel
k=k./sum(k(:));

%Convolve logic with kernel
M=double(L);
p=size(k)-1;
M_rep=zeros(size(M)+p);
M_rep(p(1)-1+(1:size(M,1)),p(2)-1+(1:size(M,2)),p(3)-1+(1:size(M,3)))=M;
Mk = convn(M_rep,k,'valid');
Mk=Mk./max(Mk(:)); %Scale max to 1

epsMax=max(eps(Mk(:)));

%Erode or dilate logic
switch erodeDilateOption
    case 'erode'
        Lk=(Mk>=(1-epsMax)); %Stayed nearly 1
    case 'dilate'    
        Lk=Mk>=(0+epsMax); %Became higher than 0
end
