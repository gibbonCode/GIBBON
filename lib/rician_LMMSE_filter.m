function [A]=rician_LMMSE_filter(M,S,k)

%Based on Aja-Fernandez et al. 2008. Restoration of DWI Data Using a Rician
%LMMSE Estimator

%%
% Sg=1; 
pad_step=(k-1)/2;
% hg=gaussiankernel(2,k,Sg); %Equivalent to: hg = fspecial('gaussian',[k k], S);

% hg=gauss_kernel(k,nd,f,m);
% hg=gauss_kernel(k,2,Sg,'sigma');
hg=gauss_kernel(k,2,1.5,'width');

Mij2=M.^2; Mij4=M.^4;

%% Gaussian filtering M.^2

[Mij2_gauss]=padrep(Mij2,pad_step.*ones(1,3)); %Padding array
Mij2_gauss=convn(Mij2_gauss,hg,'same');
Mij2_gauss=Mij2_gauss(pad_step+1:(size(Mij2_gauss,1)-pad_step),pad_step+1:(size(Mij2_gauss,2)-pad_step),pad_step+1:(size(Mij2_gauss,3)-pad_step)); %Trimming back to normal size

%% Gaussian filtering M.^4

[Mij4_gauss]=padrep(Mij4,pad_step.*ones(1,3)); %Padding array
Mij4_gauss=convn(Mij4_gauss,hg,'same');
Mij4_gauss=Mij4_gauss(pad_step+1:(size(Mij4_gauss,1)-pad_step),pad_step+1:(size(Mij4_gauss,2)-pad_step),pad_step+1:(size(Mij4_gauss,3)-pad_step)); %Trimming back to normal size

%% Calculating K

Kij_nom=(4.*(S.^2)).*(Mij2_gauss-(S.^2));
Kij_denom=Mij4_gauss-(Mij2_gauss.^2);
L=Kij_denom~=0; %Logic index of non zero values
Kij=zeros(size(Kij_denom));
Kij(L)=1-(Kij_nom(L)./Kij_denom(L));
Kij(Kij<0)=0;

%% Calculating signal A

A2= Mij2_gauss-2.*(S.^2) + Kij.*(Mij2-Mij2_gauss) ;
A2(A2<0)=0; %Smaller then zero set to zero
A=sqrt(A2);

%% END