function [D_out]=disp2strain(Ux, Uy, Uz, v, strain_type)


%% GETTING DISPLACEMENT GRADIENT MATRIX
%See help gradient; "The first output FX is always the gradient along the
%2nd dimension of F, going across columns. The second output FY is always
%the gradient along the 1st dimension of F, going across rows. For the
%third output FZ and the outputs that follow, the Nth output is the
%gradient along the Nth dimension of F."
% Therefore the gradient of A: [Gdim2,Gdim1,Gdim3] = gradient(A,v1,v2,v3)

[Ugxx,Ugxy,Ugxz] = gradient(Ux,v(1),v(2),v(3));
[Ugyx,Ugyy,Ugyz] = gradient(Uy,v(1),v(2),v(3));
[Ugzx,Ugzy,Ugzz] = gradient(Uz,v(1),v(2),v(3));

UG=[Ugxx(:) Ugxy(:) Ugxz(:)...
    Ugyx(:) Ugyy(:) Ugyz(:)...
    Ugzx(:) Ugzy(:) Ugzz(:)];
UG=reshape(UG',3,3,size(UG,1));
UG=reshape(mat2cell(UG,3,3,ones(size(UG,3),1)),size(Ux));

%% DERIVING STRAIN MEASURES
I_unity=eye(3,3);

%Allocating memory
F=UG; C=UG; E=UG; U=UG; Cp=UG;
Emaxp=ones(size(Ux)); Eminp=Emaxp; Emidp=Emaxp; Eoct=Emaxp; 
Exx=Emaxp; Exy=Emaxp; Exz=Emaxp;
Eyx=Emaxp; Eyy=Emaxp; Eyz=Emaxp;
Ezx=Emaxp; Ezy=Emaxp; Ezz=Emaxp;

%Looping to save memory
for i=1:1:size(UG,1);
    for j=1:1:size(UG,2);
        for k=1:1:size(UG,3);
            Fijk=UG{i,j,k}+I_unity;
            F(i,j,k)={Fijk};
            Cijk=Fijk'*Fijk; %The right Cauchy Green tensor
            C(i,j,k)={Cijk};
            if any(isnan(Cijk))
                Qijk=NaN(size(Cijk));
                Up_ijk=NaN(size(Cijk));
                Ep_ijk=NaN(size(Cijk));
            else
                [Qijk,Lpsq_ijk] = eig(Cijk); %The rotation (and translation) tensor and the squared principal stretches
                Up_ijk=sqrt(Lpsq_ijk); %Right stretch tensor in principal space
                Cp(i,j,k)={Lpsq_ijk};
                switch strain_type
                    case 1 %Principal Biot strain tensor
                        Ep_ijk=(Up_ijk)-I_unity; 
                    case 2 %Principal logarithmic or Hencky strain tensor                        
                        Ep_ijk=I_unity; 
                        Ep_ijk(I_unity>0)=log(diag(Up)); %Log on eigenvalues avoids log(0)=-inf related issues
                    case 3 %Principal Green-Lagrange strain tensor
                        Ep_ijk=0.5.*((Up_ijk.^2)-1); 
                end
            end
            Ep_ijk(Up_ijk==0)=0; 
            
            Uijk=Qijk*Up_ijk*Qijk'; %Right stretch tensor in cartesian space
            U(i,j,k)={Uijk};

            Eijk=Qijk*Ep_ijk*Qijk'; %Strain tensor in cartesian space
            E(i,j,k)={Eijk};

            %Principal strains
            E_diag=diag(Ep_ijk);
            [Emaxp(i,j,k),ind_max]=max(E_diag); %Maximum principal strain
            [Eminp(i,j,k),ind_min]=min(E_diag); %Minimum principal strain
            Emidp(i,j,k)=E_diag(find(~ismember([1 2 3],[ind_max ind_min]),1)); %Mid principal strain 
            %N.B. The mid principal strain is now defined as "the other
            %one". However if the min/max values are not unique these
            %functions just pick one hence mid may be equal to say min
            
            %Octahedral strain
            Eoct(i,j,k)=1/3.*sqrt(((Ep_ijk(1,1)-Ep_ijk(2,2)).^2)+((Ep_ijk(2,2)-Ep_ijk(3,3)).^2)+((Ep_ijk(3,3)-Ep_ijk(1,1)).^2)); 
            
            %Strain components
            Exx(i,j,k)=Eijk(1,1); Exy(i,j,k)=Eijk(1,2); Exz(i,j,k)=Eijk(1,3);
            Eyx(i,j,k)=Eijk(2,1); Eyy(i,j,k)=Eijk(2,2); Eyz(i,j,k)=Eijk(2,3);
            Ezx(i,j,k)=Eijk(3,1); Ezy(i,j,k)=Eijk(3,2); Ezz(i,j,k)=Eijk(3,3);
            
        end
    end
end

%Storing in structure array
D_out.F=F; 
D_out.C=C; 
D_out.Cp=Cp; 
D_out.U=U; 
D_out.E=E; 
D_out.Emaxp=Emaxp; 
D_out.Emidp=Emidp; 
D_out.Eminp=Eminp;
D_out.Eoct=Eoct;
D_out.Exx=Exx; D_out.Exy=Exy; D_out.Exz=Exz;
D_out.Eyx=Eyx; D_out.Eyy=Eyy; D_out.Eyz=Eyz;
D_out.Ezx=Ezx; D_out.Ezy=Ezy; D_out.Ezz=Ezz; 

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
