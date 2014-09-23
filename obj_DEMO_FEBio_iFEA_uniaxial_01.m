function [Fopt,OPT_stats_out]=obj_DEMO_FEBio_iFEA_uniaxial_01(Pb,objectiveStruct)

%% PARAMETER BOUNDS

P=zeros(size(Pb));
for q=1:1:numel(P);
    Psn=objectiveStruct.Ps;
    Psn.ub=objectiveStruct.Ps.ub(q);
    Psn.lb=objectiveStruct.Ps.lb(q);
    Psn.c=objectiveStruct.Ps.c(q);
    P(q)=parbound(Pb(q),Psn);
end

%% SETTING MATERIAL PARAMETERS

%Acces material parameters
mat_struct=objectiveStruct.mat_struct;
mat_struct.par_values={P(1) P(2) P(1)*objectiveStruct.k_factor}; 

disp('SETTING MATERIAL PARAMETERS...');
disp_text=sprintf('%6.16e,',P); disp_text=disp_text(1:end-1);
disp(['P=',disp_text]);

%Assign material parameters
set_mat_par_FEBIO(objectiveStruct.FEB_struct.run_filename,objectiveStruct.FEB_struct.run_filename,{mat_struct});

%% START FEBio NOW

[runFlag]=runMonitorFEBio(objectiveStruct.FEBioRunStruct);

bcPrescribeList=objectiveStruct.bcPrescribeList;
sampleHeight=objectiveStruct.sampleHeight;
initialArea=objectiveStruct.initialArea;
stretch_exp=objectiveStruct.stretch_exp;
stress_cauchy_exp=objectiveStruct.stress_cauchy_exp;

if runFlag==1  
  
    
    %Importing displacement
    [~,N_disp_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{1}); %Nodal displacements

    %Importing force
    [~,N_force_mat,~]=importFEBio_logfile(objectiveStruct.FEB_struct.run_output_names{2}); %Nodal displacements
    
    %Get Z forces
    FZ=sum(N_force_mat(bcPrescribeList,end,:),1);
    FZ=[0; FZ(:)]; %Mean top surface nodal forces
    
    %Derive applied stretch
    DZ_set=N_disp_mat(bcPrescribeList,end,:); %Final nodal displacements
    DZ_set=mean(DZ_set,1);
    stretch_sim=(DZ_set+sampleHeight)./sampleHeight;
    stretch_sim=[1; stretch_sim(:)];
    
    %Derive simulated Cauchy stress 
    currentArea=initialArea./stretch_sim;
    stress_cauchy_sim=FZ./currentArea; %Cauchy stress
    stress_cauchy_sim=stress_cauchy_sim.*1e3; %Scale to kPa
    
    %Interpolate experiment onto simulated points
    stress_cauchy_sim_exp = interp1(stretch_sim,stress_cauchy_sim,stretch_exp,'cubic');
    
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


