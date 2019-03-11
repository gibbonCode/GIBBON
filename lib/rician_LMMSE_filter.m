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
