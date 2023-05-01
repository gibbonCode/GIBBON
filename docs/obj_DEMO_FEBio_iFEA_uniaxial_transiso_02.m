function [Fopt,OPT_stats_out]=obj_DEMO_FEBio_iFEA_uniaxial_transiso_02(Pn,objectiveStruct)

%% Unnormalize and constrain parameters

P=Pn.*objectiveStruct.parNormFactors; %Scale back, undo normalization
P_in=P; %Proposed P

%Constraining parameters
for q=1:1:numel(P);
    [P(q)]=parLimNat(objectiveStruct.Pb_struct.xx_c(q),objectiveStruct.Pb_struct.xxlim(q,:),P(q));
end

%% SETTING MATERIAL PARAMETERS

mat_struct.par_names={{'solid','Ogden unconstrained','c1'},...
     {'solid','Ogden unconstrained','m1'},...
    {'solid','Ogden unconstrained','c2'},...
     {'solid','Ogden unconstrained','m2'},...
    {'solid','Ogden unconstrained','cp'},...
    {'solid','ellipsoidal fiber distribution','ksi'},...
    {'solid','ellipsoidal fiber distribution','beta'},...
    };

%Acces material parameters
mat_struct.id=1;

c1=P(1);
m1=P(2);
ksi=[P(3) P(3) P(3)*P(4)];
beta=P(5)*ones(1,3);
cp=(2*c1+mean(ksi))*objectiveStruct.k_factor;

mat_struct.par_values={c1 m1 c1 -m1 cp ksi beta};

disp('SETTING MATERIAL PARAMETERS...');
disp(['Proposed (norm.): ',sprintf(repmat('%6.16e ',[1,numel(Pn)]),Pn)]);
disp(['Proposed        : ',sprintf(repmat('%6.16e ',[1,numel(P_in)]),P_in)]);
disp(['Set (constr.)   : ',sprintf(repmat('%6.16e ',[1,numel(P)]),P)]);

%Assign material parameters
docNode=set_mat_par_FEBIO(objectiveStruct.FEB_struct.run_filename,objectiveStruct.FEB_struct.run_filename,{mat_struct});

disp('Done')

%%

for q=1:1:2 %Direction cases
    
    disp('SETTING FIBRE DIRECTIONS...');
    
    switch q
        case 1
            alphaFib=0;
        case 2
            alphaFib=0.5*pi;
    end
    
    [R,~]=euler2DCM([0,alphaFib,0]);
    v_fib=(R*[0 0 1]')';
    V_fib=v_fib(ones(numel(objectiveStruct.FEB_struct.Geometry.ElementData.MatAxis.ElementIndices),1),:);
    
    %Adding fibre direction, construct local orthonormal basis vectors
    [a,d]=vectorOrthogonalPair(V_fib);
    VF_E=zeros(size(V_fib,1),size(V_fib,2),2);
    VF_E(:,:,1)=a; %a1 ~ e1 ~ X or first direction
    VF_E(:,:,2)=d; %a2 ~ e2 ~ Y or second direction
    %Vf_E %a3 ~ e3 ~ Z, third direction, or fibre direction
    objectiveStruct.FEB_struct.Geometry.ElementData.MatAxis.Basis=VF_E;
    docNode=addMatAxisFibreElementData_FEB(docNode,objectiveStruct.FEB_struct);
    write_XML_no_extra_lines(objectiveStruct.FEB_struct.run_filename,docNode)% Saving XML file
    disp('Done')
    
    %% START FEBio NOW
    
    [runFlag]=runMonitorFEBio(objectiveStruct.FEBioRunStruct);
    
    stretch_exp=objectiveStruct.stretch_exp;
    stress_cauchy_exp=objectiveStruct.stress_cauchy_exp(:,q);
    
    if runFlag==1
        
        %Importing stress
        [~,S_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{2}); %Element Cauchy stresses
        S_mat=S_mat(:,2:end,:); %Final stress, no element labels
        S_mat=squeeze(mean(S_mat,1)); %Mean across elements
        stress_cauchy_sim=[0; S_mat(:)]; %Cauchy stress
        stress_cauchy_sim=stress_cauchy_sim.*1e3; %Scale to kPa
        
        %Import principal strains
        [~, E_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{1}); %Nodal displacements
        E_mat=E_mat(:,2:end,:);
        
        %Derive x stretch
        EX_mat=E_mat(:,1,:);
        EX_mat=squeeze(mean(EX_mat,1)); %Mean across elements
        stretch_sim_X=sqrt(2*EX_mat+1);
        stretch_sim_X=[1; stretch_sim_X(:)];
        stretch_sim_X_end=(stretch_sim_X(end));
        
        %Derive y stretch
        EY_mat=E_mat(:,2,:);
        EY_mat=squeeze(mean(EY_mat,1)); %Mean across elements
        stretch_sim_Y=sqrt(2*EY_mat+1);
        stretch_sim_Y=[1; stretch_sim_Y(:)];
        stretch_sim_Y_end=(stretch_sim_Y(end));
        
        %Derive z stretch
        EZ_mat=E_mat(:,3,:);
        EZ_mat=squeeze(mean(EZ_mat,1)); %Mean across elements
        stretch_sim_Z=sqrt(2*EZ_mat+1);
        stretch_sim_Z=[1; stretch_sim_Z(:)];
        stretch_sim_Z_end=(stretch_sim_Z(end));
        
        %Interpolate experiment onto simulated points
        stress_cauchy_sim_exp = interp1(stretch_sim_Z,stress_cauchy_sim,stretch_exp,'pchip');
        
        %Derive Fopt
        stressDev=stress_cauchy_exp-stress_cauchy_sim_exp;        
        Fopt_stress=mean(abs(stressDev)./max(abs(stress_cauchy_exp)));

        switch q
            case 1
                stretchDev=[1./(sqrt(0.7)) 1./(sqrt(0.7))]-[stretch_sim_X_end stretch_sim_Y_end];
            case 2
                stretchDev=[exp(log(stretch_sim_Z_end)*-0.36) exp(log(stretch_sim_Z_end)*-0.65)]-[stretch_sim_X_end stretch_sim_Y_end];
        end
        
        Fopt_stretch=mean(abs(stretchDev));
        
        [R_sq]=R_squared(stress_cauchy_sim_exp,stress_cauchy_exp);

    else %Output NaN
        stress_cauchy_sim=NaN(size(stretch_exp));
        stretch_sim_Z=NaN(size(stretch_exp));
        stressDev=NaN(size(stretch_exp));
        stretchDev=NaN(1,2);
        Fopt_stress=Nan;  
        Fopt_stretch=NaN;
    end
    
    OPT_stats_out.stress_cauchy_sim{q}=stress_cauchy_sim;
    OPT_stats_out.stretch_sim{q}=stretch_sim_Z;
    OPT_stats_out.stretch_sim_end(q,:)=[stretch_sim_X_end stretch_sim_Y_end stretch_sim_Z_end];
    OPT_stats_out.stressDev{q}=stressDev;
    OPT_stats_out.stretchDev{q}=stretchDev;
    OPT_stats_out.Fopt(:,q)=[Fopt_stress Fopt_stretch];
    OPT_stats_out.R_sq(q)=R_squared(stress_cauchy_exp,stress_cauchy_sim_exp);
    
end

switch objectiveStruct.method
    case 1
        Fopt=OPT_stats_out.Fopt(1,1)+OPT_stats_out.Fopt(1,2)+OPT_stats_out.Fopt(2,2); %Scalar objective function
    case 2
        Fopt=[OPT_stats_out.Fopt(1,1) OPT_stats_out.Fopt(1,2) OPT_stats_out.Fopt(2,2)]; %Objective function vectors
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
