function [Fopt,OPT_stats_out]=obj_DEMO_FEBio_iFEA_uniaxial_transiso_01(Pn,objectiveStruct)

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
        {'solid','fiber-exp-pow','ksi'},...
        {'solid','fiber-exp-pow','beta'},...
        };
    
%Acces material parameters
mat_struct.id=1;
if objectiveStruct.dirCase==1
    alphaFib=0;   
    c1=P(1);
    m1=P(2);
    ksi=objectiveStruct.P_ini(3);
    beta=objectiveStruct.P_ini(4);
    cp=(c1+ksi)*objectiveStruct.k_factor;
elseif objectiveStruct.dirCase==2
    alphaFib=0.5*pi;
    c1=objectiveStruct.P_ini(1); 
    m1=objectiveStruct.P_ini(2);
    ksi=P(1);
    beta=P(2);
    cp=(c1+ksi)*objectiveStruct.k_factor;    
end
mat_struct.par_values={c1 m1 c1 -m1 cp ksi beta};

disp('SETTING MATERIAL PARAMETERS...');
disp(['Proposed (norm.): ',sprintf(repmat('%6.16e ',[1,numel(Pn)]),Pn)]);
disp(['Proposed        : ',sprintf(repmat('%6.16e ',[1,numel(P_in)]),P_in)]);
disp(['Set (constr.)   : ',sprintf(repmat('%6.16e ',[1,numel(P)]),P)]);

%Assign material parameters
docNode=set_mat_par_FEBIO(objectiveStruct.FEB_struct.run_filename,objectiveStruct.FEB_struct.run_filename,{mat_struct});

disp('Done')

%%

disp('SETTING FIBRE DIRECTIONS...');

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
docNode=add_fiber_dir_FEB(docNode,objectiveStruct.FEB_struct,0);
write_XML_no_extra_lines(objectiveStruct.FEB_struct.run_filename,docNode)% Saving XML file
disp('Done')

%% START FEBio NOW

[runFlag]=runMonitorFEBio(objectiveStruct.FEBioRunStruct);

stretch_exp=objectiveStruct.stretch_exp;
stress_cauchy_exp=objectiveStruct.stress_cauchy_exp;

if runFlag==1
    
    %Importing stress
    [~,S_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{2}); %Element Cauchy stresses
    S_mat=S_mat(:,2:end,:); %Final stress, no element labels
    S_mat=squeeze(mean(S_mat,1)); %Mean across elements
    stress_cauchy_sim=[0; S_mat(:)]; %Cauchy stress
    stress_cauchy_sim=stress_cauchy_sim.*1e3; %Scale to kPa
    
    %Import displacement
    [~, N_disp_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{1}); %Nodal displacements
    
    %Derive applied stretch
    DZ_set=N_disp_mat(objectiveStruct.bcPrescribeList,end,:); %Final nodal displacements
    DZ_set=mean(DZ_set,1);
    stretch_sim=(DZ_set+objectiveStruct.sampleHeight)./objectiveStruct.sampleHeight;
    stretch_sim=[1; stretch_sim(:)];
    
    %Interpolate experiment onto simulated points
    stress_cauchy_sim_exp = interp1(stretch_sim,stress_cauchy_sim,stretch_exp,'pchip');
    
    %Derive Fopt
    stressDev=stress_cauchy_exp-stress_cauchy_sim_exp;
    
    switch objectiveStruct.method
        case 1
            Fopt=sum((stressDev).^2); %Sum of squared differences
        case 2
            Fopt=(stressDev).^2; %Squared differences
    end
    
    OPT_stats_out.stress_cauchy_sim=stress_cauchy_sim;
    OPT_stats_out.stretch_sim=stretch_sim;
    OPT_stats_out.stressDev=stressDev;
    OPT_stats_out.Fopt=Fopt;
    
else %Output NaN
    switch objectiveStruct.method
        case 1
            Fopt=NaN;
        case 2
            Fopt=NaN(size(stress_cauchy_exp)); %Squared differences
    end
    OPT_stats_out=[];
end

%%


