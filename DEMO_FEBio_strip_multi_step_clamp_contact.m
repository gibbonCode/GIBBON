%% DEMO_FEBio_strip_multi_step_clamp_contact
% Below is a demonstration for: 
% 1) The creation of an FEBio model for clamped tensile testing
% 2) The use of multiple steps
% 4) Running an FEBio job with MATLAB
% 5) Importing FEBio results into MATLAB

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=50;

%%
% Control parameters

% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath),'data','temp');

modelName='tempModel';
modelNameFull=fullfile(savePath,modelName);

pointspacing=0.75; 

%Specifying dimensions and number of elements
sampleWidth=10;
numElementsWidth=round(sampleWidth/pointspacing);
numElementsWidth=numElementsWidth+iseven(numElementsWidth); %Force uneven so there is a middle element
elementSizeWidth=sampleWidth/numElementsWidth;

sampleThickness=1.65; 
numElementsThickness=round(sampleThickness/pointspacing)+1;

sampleGripGripHeight=20;
numElementsGripGripHeight=round(sampleGripGripHeight/pointspacing);
numElementsGripGripHeight=numElementsGripGripHeight+iseven(numElementsGripGripHeight); %Force uneven so there is a middle element

sampleClampedHeight=8;
numElementsClampedHeight=round(sampleClampedHeight/pointspacing);
elementSizeClamped=sampleClampedHeight/numElementsClampedHeight;

contactOverlap=3;

contactInitialOffset=0.1;
clampCompressiveStrain=0.375;
clampCompressiveDisplacement=(sampleThickness.*clampCompressiveStrain)/2;
tensileStretch=1.1;
clampTensionDisplacement=(sampleGripGripHeight.*tensileStretch)-sampleGripGripHeight;

contactType=1; %1=sticky, 2=sliding/friction

%% CREATING 3 MESHED BOXES

%Create box 1
boxDim=[sampleWidth sampleThickness sampleClampedHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsClampedHeight]; %Number of elements
[box1]=hexMeshBox(boxDim,boxEl);
E1=box1.E;
V1=box1.V;
F1=box1.F;
Fb1=box1.Fb;
faceBoundaryMarker1=box1.faceBoundaryMarker;

%Create box 3 by copying the first
E3=E1; 
V3=V1; 
F3=F1;
Fb3=Fb1;
faceBoundaryMarker3=faceBoundaryMarker1;

%Shift first box up
V1(:,3)=V1(:,3)+sampleGripGripHeight/2+sampleClampedHeight/2;

%Shift third box down
V3(:,3)=V3(:,3)-sampleGripGripHeight/2-sampleClampedHeight/2;

%Create box 1
boxDim=[sampleWidth sampleThickness sampleGripGripHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsGripGripHeight]; %Number of elements
[box2]=hexMeshBox(boxDim,boxEl);
E2=box2.E;
V2=box2.V;
F2=box2.F;
Fb2=box2.Fb;
faceBoundaryMarker2=box2.faceBoundaryMarker;

%%
% Plotting surface models
hf=cFigure;
title('Box sets','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb1,'Vertices',V1,'FaceColor','flat','CData',faceBoundaryMarker1,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Fb2,'Vertices',V2,'FaceColor','flat','CData',faceBoundaryMarker2,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Fb3,'Vertices',V3,'FaceColor','flat','CData',faceBoundaryMarker3,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
camlight headlight;
drawnow; 

%% MERGING BOX SETS

faceBoundaryMarker_all=[faceBoundaryMarker1; faceBoundaryMarker2; faceBoundaryMarker3;];
faceBoundaryMarker_ind=[ones(size(Fb1,1),1);2*ones(size(Fb2,1),1); 3*ones(size(Fb3,1),1);];

V=[V1;V2;V3];
E=[E1;E2+size(V1,1);E3+size(V1,1)+size(V2,1)];
F=[F1;F2+size(V1,1);F3+size(V1,1)+size(V2,1)];
Fb=[Fb1;Fb2+size(V1,1);Fb3+size(V1,1)+size(V2,1)];

[~,ind1,ind2]=unique(pround(V,5),'rows');
V=V(ind1,:);
E=ind2(E);
F=ind2(F);
Fb=ind2(Fb);

%%
% Plotting surface models
hf=cFigure;
title('Merged box sets','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker_all,'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
drawnow; 

%% Create grip models

%Get edge lengths of model
% [D]=patchEdgeLengths(Fb,V);
% elementSizeRigid=min(D)/2;

elementSizeRigid=min(elementSizeClamped,elementSizeWidth)/2;

nX=round((2*sampleWidth+2*contactOverlap)/elementSizeRigid);
xRange=linspace(-(sampleWidth/2)-contactOverlap,(sampleWidth/2)+contactOverlap,nX);

nZ=round((sampleClampedHeight+contactOverlap)/elementSizeRigid);
zRange=linspace(0,sampleClampedHeight+contactOverlap,nZ);

[Xr,Zr]=meshgrid(xRange,zRange); %mesh of single slice
[Fr,Vr,~] = surf2patch(Xr,zeros(size(Xr)),Zr,zeros(size(Xr))); %Convert to patch data (quadrilateral faces)

%Create off set rigid parts
Fr1=fliplr(Fr);
Vr1=Vr;
Vr1(:,3)=Vr1(:,3)+sampleGripGripHeight/2;
Vr1(:,2)=Vr1(:,2)-sampleThickness/2-contactInitialOffset;

Fr2=Fr;
Vr2=Vr;
Vr2(:,3)=Vr2(:,3)+sampleGripGripHeight/2;
Vr2(:,2)=Vr2(:,2)+sampleThickness/2+contactInitialOffset;

Fr3=Fr;
Vr3=-Vr1;

Fr4=fliplr(Fr);
Vr4=-Vr2;

%%
% Plotting surface models
hf=cFigure;
title('Sample with rigid body clamps','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',F,'Vertices',V,'FaceColor',0.5*ones(1,3),'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

patch('Faces',Fr1,'Vertices',Vr1,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fr1,Vr1,2);
patch('Faces',Fr2,'Vertices',Vr2,'FaceColor','g','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fr2,Vr2,2);
patch('Faces',Fr3,'Vertices',Vr3,'FaceColor','b','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fr3,Vr3,2);
patch('Faces',Fr4,'Vertices',Vr4,'FaceColor','y','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fr4,Vr4,2);
set(gca,'FontSize',fontSize);

view(3); axis tight;  axis equal;  grid on;
drawnow; 

%% Define contact surfaces

logicContactSurf1=faceBoundaryMarker_all==3 & faceBoundaryMarker_ind==1;
Fc1=Fb(logicContactSurf1,:);

logicContactSurf2=faceBoundaryMarker_all==4 & faceBoundaryMarker_ind==1;
Fc2=Fb(logicContactSurf2,:);

logicContactSurf3=faceBoundaryMarker_all==4 & faceBoundaryMarker_ind==3;
Fc3=Fb(logicContactSurf3,:);

logicContactSurf4=faceBoundaryMarker_all==3 & faceBoundaryMarker_ind==3;
Fc4=Fb(logicContactSurf4,:);

% Plotting surface models
hf=cFigure;
title('Contact sets','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fc1,'Vertices',V,'FaceColor','r','FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fc1,V,2);
patch('Faces',Fc2,'Vertices',V,'FaceColor','g','FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fc2,V,2);
patch('Faces',Fc3,'Vertices',V,'FaceColor','b','FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fc3,V,2);
patch('Faces',Fc4,'Vertices',V,'FaceColor','y','FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(Fc4,V,2);

patch('Faces',Fr1,'Vertices',Vr1,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Fr2,'Vertices',Vr2,'FaceColor','g','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Fr3,'Vertices',Vr3,'FaceColor','b','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Fr4,'Vertices',Vr4,'FaceColor','y','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);

view(3); axis tight;  axis equal;  grid on;
drawnow; 


%% Collect node sets

VT=[V; Vr1; Vr2; Vr3; Vr4];

%Fix rigid model node indices
Eq1=Fr1+size(V,1);
Eq2=Fr2+size(V,1)+size(Vr1,1);
Eq3=Fr3+size(V,1)+size(Vr1,1)+size(Vr2,1);
Eq4=Fr4+size(V,1)+size(Vr1,1)+size(Vr2,1)+size(Vr3,1);

bcPrescribeListRB1=unique(Eq1(:));
bcPrescribeMagnitudesRB1=zeros(1,3);
bcPrescribeMagnitudesRB1(:,1)=0;
bcPrescribeMagnitudesRB1(:,2)=clampCompressiveDisplacement+contactInitialOffset; %In step 1
bcPrescribeMagnitudesRB1(:,3)=clampTensionDisplacement; %In step 2

bcPrescribeListRB2=unique(Eq2(:));
bcPrescribeMagnitudesRB2=zeros(1,3);
bcPrescribeMagnitudesRB2(:,1)=0;
bcPrescribeMagnitudesRB2(:,2)=-clampCompressiveDisplacement-contactInitialOffset; %In step 1
bcPrescribeMagnitudesRB2(:,3)=clampTensionDisplacement; %In step 2

bcPrescribeListRB3=unique(Eq3(:));
bcPrescribeMagnitudesRB3=zeros(1,3);
bcPrescribeMagnitudesRB3(:,1)=0;
bcPrescribeMagnitudesRB3(:,2)=-clampCompressiveDisplacement-contactInitialOffset; %In step 1
bcPrescribeMagnitudesRB3(:,3)=0; %In step 2

bcPrescribeListRB4=unique(Eq4(:));
bcPrescribeMagnitudesRB4=zeros(1,3);
bcPrescribeMagnitudesRB4(:,1)=0;
bcPrescribeMagnitudesRB4(:,2)=clampCompressiveDisplacement+contactInitialOffset; %In step 1
bcPrescribeMagnitudesRB4(:,3)=0; %In step 2

% Plotting surface models
hf=cFigure;
title('Complete model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',VT,'FaceColor',0.5*ones(1,3),'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

patch('Faces',Eq1,'Vertices',VT,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Eq2,'Vertices',VT,'FaceColor','g','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Eq3,'Vertices',VT,'FaceColor','b','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
patch('Faces',Eq4,'Vertices',VT,'FaceColor','y','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);

set(gca,'FontSize',fontSize);

view(3); axis tight;  axis equal;  grid on;
drawnow; 

%% CONSTRUCTING FEB MODEL

% Defining file names
FEB_struct.run_filename=[modelNameFull,'.feb']; %FEB file name
FEB_struct.run_logname=[modelNameFull,'.txt']; %FEBio log file name

%Creating FEB_struct
FEB_struct.Geometry.Nodes=VT;
FEB_struct.Geometry.Elements={Eq1 Eq2 Eq3 Eq4 E}; %The element sets
FEB_struct.Geometry.ElementType={'quad4','quad4','quad4','quad4','hex8'}; %The element types
FEB_struct.Geometry.ElementMat={[1*ones(1,size(Eq1,1))]; [2*ones(1,size(Eq2,1))]; [3*ones(1,size(Eq3,1))]; [4*ones(1,size(Eq4,1))]; [5*ones(1,size(E,1))];};

%Adding fibre direction, construct local orthonormal basis vectors
Vf_E=zeros(size(E,1),3);
Vf_E(:,3)=1; %Z-axis fibres
% Vf_E(:,1)=1; %X-axis fibres
[a,d]=vectorOrthogonalPair(Vf_E);

VF_E=nan(size(Vf_E,1),size(Vf_E,2),2);
VF_E(:,:,1)=a; %a1 ~ e1 ~ X or first direction
VF_E(:,:,2)=d; %a2 ~ e2 ~ Y or second direction
%Vf_E %a3 ~ e3 ~ Z, third direction, or fibre direction

FEB_struct.Geometry.ElementData.MatAxis.ElementIndices=size(Eq1,1)+size(Eq2,1)+size(Eq3,1)+size(Eq4,1)+(1:1:size(E,1));
FEB_struct.Geometry.ElementData.MatAxis.Basis=VF_E;

% DEFINING MATERIALS

%Rigid materials
Mat1.type='rigid body';
Mat1.props={'density','center_of_mass'};
Mat1.vals={1,[0,0,0]};
Mat1.aniso_type='none';

Mat2.type='rigid body';
Mat2.props={'density','center_of_mass'};
Mat2.vals={1,[0,0,0]};
Mat2.aniso_type='none';

Mat3.type='rigid body';
Mat3.props={'density','center_of_mass'};
Mat3.vals={1,[0,0,0]};
Mat3.aniso_type='none';

Mat4.type='rigid body';
Mat4.props={'density','center_of_mass'};
Mat4.vals={1,[0,0,0]};
Mat4.aniso_type='none';

% c1=1e-3;
% k=c1*100;
% Mat5.type='Mooney-Rivlin';
% Mat5.props={'c1','c2','k'};
% Mat5.vals={c1,0,k};
% Mat5.aniso_type='none';

c1=2.309;
m1=9.421;
ksi=22.499;
beta=2.387;
k=(0.5.*(c1+ksi))*100;% 20;
Mat5.type='uncoupled solid mixture';
Mat51.type='Ogden';
Mat51.props={'c1','m1','k'};
Mat51.vals={c1,m1,k};
Mat51.aniso_type='none';
Mat52.type='fiber-exp-pow-uncoupled';
Mat52.props={'ksi','alpha','beta','theta','phi','k'};
Mat52.vals={ksi,1e-25,beta,0,0,k};
Mat52.aniso_type='none';
Mat5.Mats={Mat51 Mat52};

FEB_struct.Materials={Mat1 Mat2 Mat3 Mat4 Mat5};


%Step specific BC's
% FEB_struct.Step(1).Boundary.PrescribeList={bcPrescribeList,bcPrescribeList,bcPrescribeList};
% FEB_struct.Step(1).Boundary.PrescribeType={'x','y','z'};
% FEB_struct.Step(1).Boundary.PrescribeValues={bcPrescribedMagnitudesStep1(:,1),bcPrescribedMagnitudesStep1(:,2),bcPrescribedMagnitudesStep1(:,3)};
% FEB_struct.Step(1).Boundary.PrescribeTypes={'relative','relative','relative'};
% FEB_struct.Step(1).Boundary.LoadCurveIds=[1 1 1];
% 
% FEB_struct.Step(2).Boundary.PrescribeList={bcPrescribeList,bcPrescribeList,bcPrescribeList};
% FEB_struct.Step(2).Boundary.PrescribeType={'x','y','z'};
% FEB_struct.Step(2).Boundary.PrescribeValues={bcPrescribedMagnitudesStep2(:,1),bcPrescribedMagnitudesStep2(:,2),bcPrescribedMagnitudesStep2(:,3)};
% FEB_struct.Step(2).Boundary.PrescribeTypes={'relative','relative','relative'};
% FEB_struct.Step(2).Boundary.LoadCurveIds=[2 2 2];

%Step specific control sections
FEB_struct.Step(1).Control.AnalysisType='static';
FEB_struct.Step(1).Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Step(1).Control.Values={10,0.1,...
    25,0,...
    0.001,0.01,0,0.9};
FEB_struct.Step(1).Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Step(1).Control.TimeStepperValues={1e-5, 0.1, 5, 5, 1};

FEB_struct.Step(2).Control=FEB_struct.Step(1).Control;

%Constraint section
FEB_struct.Constraints(1).RigidId=1;
FEB_struct.Constraints(1).Properties={'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
FEB_struct.Constraints(1).Values={bcPrescribeMagnitudesRB1(1),bcPrescribeMagnitudesRB1(2),bcPrescribeMagnitudesRB1(3),0,0,0};
FEB_struct.Constraints(1).LoadCurveIds=[1 1 2 1 1 1];
FEB_struct.Constraints(1).Type='prescribed';

FEB_struct.Constraints(2).RigidId=2;
FEB_struct.Constraints(2).Properties={'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
FEB_struct.Constraints(2).Values={bcPrescribeMagnitudesRB2(1),bcPrescribeMagnitudesRB2(2),bcPrescribeMagnitudesRB2(3),0,0,0};
FEB_struct.Constraints(2).LoadCurveIds=[1 1 2 1 1 1];
FEB_struct.Constraints(2).Type='prescribed';

FEB_struct.Constraints(3).RigidId=3;
FEB_struct.Constraints(3).Properties={'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
FEB_struct.Constraints(3).Values={bcPrescribeMagnitudesRB3(1),bcPrescribeMagnitudesRB3(2),bcPrescribeMagnitudesRB3(3),0,0,0};
FEB_struct.Constraints(3).LoadCurveIds=[1 1 2 1 1 1];
FEB_struct.Constraints(3).Type='prescribed';

FEB_struct.Constraints(4).RigidId=4;
FEB_struct.Constraints(4).Properties={'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
FEB_struct.Constraints(4).Values={bcPrescribeMagnitudesRB4(1),bcPrescribeMagnitudesRB4(2),bcPrescribeMagnitudesRB4(3),0,0,0};
FEB_struct.Constraints(4).LoadCurveIds=[1 1 2 1 1 1];
FEB_struct.Constraints(4).Type='prescribed';

%Adding contact information
FEB_struct.Boundary.Contact{1}.Surfaces.elements={Eq1,Fc1};
FEB_struct.Boundary.Contact{1}.Surfaces.Type={'master','slave'};
switch contactType
    case 1 %STICKY
        %           <contact type="sticky">
        % 			<laugon>0</laugon>
        % 			<tolerance>0.2</tolerance>
        % 			<penalty>1</penalty>
        % 			<minaug>0</minaug>
        % 			<maxaug>10</maxaug>
        
        FEB_struct.Boundary.Contact{1}.Type='sticky';
        FEB_struct.Boundary.Contact{1}.Properties={'laugon','tolerance','penalty',...
            'minaug','maxaug'};
        FEB_struct.Boundary.Contact{1}.Values={0,0.05,(0.5.*(c1+ksi))*100,...
            0,10};
    case 2 %SLIDING/FRICTION
        FEB_struct.Boundary.Contact{1}.Type='sliding_with_gaps';
        FEB_struct.Boundary.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'fric_coeff','fric_penalty',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Boundary.Contact{1}.Values={50,0,0,...
            0,0.1,...
            0,0,10,...
            1,1000,...
            0,...
            0.01};
end

FEB_struct.Boundary.Contact{2}=FEB_struct.Boundary.Contact{1};
FEB_struct.Boundary.Contact{2}.Surfaces.elements={Eq2,Fc2};

FEB_struct.Boundary.Contact{3}=FEB_struct.Boundary.Contact{1};
FEB_struct.Boundary.Contact{3}.Surfaces.elements={Eq3,Fc3};

FEB_struct.Boundary.Contact{4}=FEB_struct.Boundary.Contact{1};
FEB_struct.Boundary.Contact{4}.Surfaces.elements={Eq4,Fc4};

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness','contact force','reaction forces'};

%Specify log file output
run_output_name_disp=[modelName,'_node_out.txt'];
run_output_name_force=[modelName,'_force_out.txt'];
FEB_struct.run_output_names={run_output_name_disp,run_output_name_force};
FEB_struct.output_types={'node_data','rigid_body_data'};
FEB_struct.data_types={'ux;uy;uz','Fx;Fy;Fz'};

%Load curves
FEB_struct.LoadData.LoadCurves.id=[1 2];
FEB_struct.LoadData.LoadCurves.type={'linear','linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1];[0 0;1 0;2 1];};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars option
febStruct2febFile_v1p2(FEB_struct);

%% RUNNING FEBIO JOB

FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1;
FEBioRunStruct.disp_log_on=1;
FEBioRunStruct.runMode='external';%'internal';
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!

%% IMPORTING NODAL DISPLACEMENT RESULTS
% Importing nodal displacements from a log file
[~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements

DN=N_disp_mat(:,2:end,end); %Final nodal displacements

%% CREATING NODE SET IN DEFORMED STATE
VT_def=VT+DN;
DN_magnitude=sqrt(sum(DN.^2,2));

%%
% Plotting the deformed model

[CF]=vertexToFaceMeasure(F,DN_magnitude);

hf1=cFigure;
title('The deformed model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',F,'Vertices',VT_def,'FaceColor','flat','CData',CF);

view(3); axis tight;  axis equal;  grid on;
colormap jet; colorbar;
% camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
