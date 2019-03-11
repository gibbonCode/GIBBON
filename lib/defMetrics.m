function [D_out]=defMetrics(F_cell,strain_type)

siz=size(F_cell);

if length(siz)==2
    siz(3)=1;
end

%% Allocating memory

%Scalar quantaties
Emaxp=ones(siz);
Eminp=Emaxp;
Emidp=Emaxp;

Lmaxp=Emaxp;
Lminp=Emaxp;
Lmidp=Emaxp;

Eoct=Emaxp;

Exx=Emaxp; Exy=Emaxp; Exz=Emaxp;
Eyx=Emaxp; Eyy=Emaxp; Eyz=Emaxp;
Ezx=Emaxp; Ezy=Emaxp; Ezz=Emaxp;

%Tensor quantaties
C=F_cell;
Cp=F_cell;
U=F_cell;
E=F_cell;
Q=F_cell;
Qeig=F_cell;

%% Looping to save memory

%Looping to save memory
for q=1:1:numel(F_cell);
    
    Fijk=F_cell{q}; %The deformation gradient tensor
    Cijk=Fijk'*Fijk; %The right Cauchy Green tensor
    C{q}=Cijk;
    if any(isnan(Cijk(:)))
        Qijk=NaN(size(Cijk));
        Up_ijk=NaN(size(Cijk));
        Ep_ijk=NaN(size(Cijk));
    else
        [Qijk,Lpsq_ijk] = eig(Cijk); %The rotation (and translation) tensor and the squared principal stretches                
        Up_ijk=sqrt(Lpsq_ijk); %Right stretch tensor in principal space
        Cp{q}=Lpsq_ijk;
        switch strain_type
            case 1 %Principal infinitesimal strain tensor
                Ep_ijk=(Up_ijk)-1;
            case 2 %Principal logarithmic strain tensor
                Ep_ijk=log(Up_ijk);
            case 3 %Principal Green-Lagrange strain tensor
                Ep_ijk=0.5.*((Up_ijk.^2)-1);
        end
    end
    Ep_ijk(Up_ijk==0)=0;
    
    Uijk=Qijk*Up_ijk*Qijk'; %Right stretch tensor in cartesian space
    U{q}=Uijk;
    
    Q{q}=Qijk;
    Qeig_ijk=Fijk/Uijk;
    Qeig{q}=Qeig_ijk;
    
    Eijk=Qijk*Ep_ijk*Qijk'; %Strain tensor in cartesian space
    E{q}=Eijk;
    
    %Principal strains
    E_diag=diag(Ep_ijk);
    [Emaxp(q),ind_max]=max(E_diag); %Maximum principal strain
    [Eminp(q),ind_min]=min(E_diag); %Minimum principal strain
    Emidp(q)=E_diag(find(~ismember([1 2 3],[ind_max ind_min]),1)); %Mid principal strain
    %N.B. The mid principal strain is now defined as "the other
    %one". However if the min/max values are not unique these
    %functions just pick one hence mid may be equal to say min
    
    %Principal stretches
    L_diag=diag(Up_ijk);
    [Lmaxp(q),ind_max]=max(L_diag); %Maximum principal stretch
    [Lminp(q),ind_min]=min(L_diag); %Minimum principal stretch
    Lmidp(q)=L_diag(find(~ismember([1 2 3],[ind_max ind_min]),1)); %Mid principal stretch
    
    %Octahedral strain
    Eoct(q)=1/3.*sqrt(((Ep_ijk(1,1)-Ep_ijk(2,2)).^2)+((Ep_ijk(2,2)-Ep_ijk(3,3)).^2)+((Ep_ijk(3,3)-Ep_ijk(1,1)).^2));
    
    %Strain components
    Exx(q)=Eijk(1,1); Exy(q)=Eijk(1,2); Exz(q)=Eijk(1,3);
    Eyx(q)=Eijk(2,1); Eyy(q)=Eijk(2,2); Eyz(q)=Eijk(2,3);
    Ezx(q)=Eijk(3,1); Ezy(q)=Eijk(3,2); Ezz(q)=Eijk(3,3);
    
end

%% Storing in structure array
D_out.F=F_cell;
D_out.C=C;
D_out.Cp=Cp;
D_out.U=U;
D_out.E=E;
D_out.Q=Q;
D_out.Qeig=Qeig;
D_out.Emaxp=Emaxp;
D_out.Emidp=Emidp;
D_out.Eminp=Eminp;
D_out.Lmaxp=Lmaxp;
D_out.Lmidp=Lmidp;
D_out.Lminp=Lminp;
D_out.Eoct=Eoct;
D_out.Exx=Exx; D_out.Exy=Exy; D_out.Exz=Exz;
D_out.Eyx=Eyx; D_out.Eyy=Eyy; D_out.Eyz=Eyz;
D_out.Ezx=Ezx; D_out.Ezy=Ezy; D_out.Ezz=Ezz;
 
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
