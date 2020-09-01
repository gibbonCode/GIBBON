%% DEMO_febio_0036_sphere_branches_slide_body_force
% Below is a demonstration for:
% 
% * Building geometry for a spherical blob with hexahedral elements
% which is pushed through a branched network. 
% This demo consists off:
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * indentation
% * contact, sliding, friction
% * rigid body constraints
% * tetrahedral elements, tet4
% * triangular elements, tri3
% * shell elements
% * sphere
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
markerSize=40;
lineWidth=3;
cMap=blood(250);

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
% febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

numElementsMantel=5;
sphereRadius=3;
voxelSize=(sphereRadius/3.5)*ones(1,3); 
w=sphereRadius*1.5;
reduceFactor=2.75;
contourLevel=1; 

% Material parameter set
c1=1e-3; %Shear-modulus-like parameter MPa
m1=2; %Material parameter setting degree of non-linearity
k_factor=10; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
d=1e-9; %Density
viscoElastic=1;
if viscoElastic==1    
    g1=2; %Viscoelastic QLV proportional coefficient
    t1=150; %Viscoelastic QLV time coefficient
end

% FEA control settings
analysisType='dynamic';
timeTotal=2; %Analysis time
numTimeSteps=200; %Number of time steps desired
step_size=timeTotal/numTimeSteps;
dtmin=(timeTotal/numTimeSteps)/100; %Minimum time step size
dtmax=timeTotal/50; %Maximum time step size
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=25; %Maximum number of retires
symmetric_stiffness=0;
min_residual=1e-20;

%Contact parameters
contactPenalty=20;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.25; 

%Specifying load
sphereVolume=4/3*(pi*sphereRadius^3); %Sphere Volume in mm^3
sphereMass=sphereVolume.*d; %Sphere mass in tone
sphereSectionArea=pi*sphereRadius^2;
bodyLoadMagnitude=(9.81*1000)*10;

forceBodyLoad=sphereMass.*bodyLoadMagnitude;
stressBodyLoad=forceBodyLoad/sphereSectionArea;

%% Define branch curves

n=100;

V1=[0 0 0; -2*sphereRadius 0 0];
V1=evenlySampleCurve(V1,n,'linear',0);

x=linspace(1,-1,n)';
xx=nan(size(x));
f=4;
if f==0
    xx=x;
else
    xx(x>=0) =(exp(-f*abs(x(x>0)))-1)./(exp(-f)-1);
    xx(x<0)  =-(exp(-f*abs(x(x<0)))-1)./(exp(-f)-1);
end
xx=w/2.*xx;
xx=-xx+w/2;

x=x.*w;
x=x+(min(V1(:,1)))+(-w);

V2=[x zeros(size(x)) xx];
V2=evenlySampleCurve(V2,n,'linear',0);

x=linspace(2,-2,n)';
xx=nan(size(x));
f=4;
if f==0
    xx=x;
else
    xx(x>=0) =(exp(-f*abs(x(x>0)))-1)./(exp(-f)-1);
    xx(x<0)  =-(exp(-f*abs(x(x<0)))-1)./(exp(-f)-1);
end
xx=w/4.*xx;
xx=-xx+w/4;

x=x.*w/2;
x=x+(-w/2);

V3=[x zeros(size(x)) xx];
V3=V3-V3(1,:);
V4=V3;
V3=V3+V2(end,:);
V3=evenlySampleCurve(V3,n,'linear',0);

V4(:,3)=-V4(:,3)/3;
V4=V4+V2(end,:);
V4=evenlySampleCurve(V4,n,'linear',0);

V5=V2;
V5(:,3)=-V5(:,3)/3;

V6=V4;
V6=V6-V6(1,:);
V6(:,3)=-V6(:,3);
V6=V6+V5(end,:);

V7=V3;
V7=V7-V7(1,:);
V7(:,3)=-V7(:,3);
V7=V7+V5(end,:);

%%

cFigure; hold on;
plotV(V1,'r.-','MarkerSize',15);
plotV(V2,'g.-','MarkerSize',15);
plotV(V3,'y.-','MarkerSize',15);
plotV(V4,'k.-','MarkerSize',15);

plotV(V5,'b.-','MarkerSize',15);
plotV(V6,'c.-','MarkerSize',15);
plotV(V7,'m.-','MarkerSize',15);
axisGeom; 
drawnow;

%%

r1=sphereRadius; %Inlet radius
A1=r1^2*pi; %Inlet area

A2=(A1*0.8); %First bifurcation inlet area
r2=sqrt(A2/pi); %First bifurcation inlet radius

A3=(A2/2); %area
r3=sqrt(A3/pi); %radius

A4=(A3/3); %area
r4=sqrt(A4/pi); %radius

E1=[(1:size(V1,1)-1)' (2:size(V1,1))'];
E2=[(1:size(V2,1)-1)' (2:size(V2,1))'];
E3=[(1:size(V3,1)-1)' (2:size(V3,1))'];
E4=[(1:size(V4,1)-1)' (2:size(V4,1))'];
E5=[(1:size(V5,1)-1)' (2:size(V5,1))'];
E6=[(1:size(V6,1)-1)' (2:size(V6,1))'];
E7=[(1:size(V7,1)-1)' (2:size(V7,1))'];

C=[linspace(r1,r2,size(V1,1))';...
   linspace(r2,r3,size(V2,1))';...
   linspace(r3,r4,size(V3,1))';...
   linspace(r3,r4,size(V4,1))';...
   linspace(r2,r3,size(V5,1))';...
   linspace(r3,r4,size(V6,1))';...
   linspace(r3,r4,size(V7,1))'];

[E,V]=joinElementSets({E1,E2,E3,E4,E5,E6,E7},{V1,V2,V3,V4,V5,V6,V7});

[E,V,ind1]=mergeVertices(E,V);
C=C(ind1);

% cPar.n=25;
% [C]=patchSmooth(E,C,[],cPar);

%%
 
cFigure; 
gpatch(E,V,'none',C,0,3);
axisGeom;
drawnow; 

%%

x=min(V(:,1)):voxelSize:max(V(:,1));
y=min(V(:,2))-1.5*sphereRadius:voxelSize:max(V(:,2))+1.5*sphereRadius;
z=min(V(:,3))-1.5*sphereRadius:voxelSize:max(V(:,3))+1.5*sphereRadius;
[X,Y,Z]=ndgrid(x,y,z);
V_grid=[X(:) Y(:) Z(:)];

[M_grid,indMin]=minDist(V_grid,V);

M_grid=M_grid./C(indMin);

M=reshape(M_grid,size(X));

%%

vizStruct.colormap=warmcold(250);
vizStruct.clim=[0 contourLevel*2];
sv3(M,voxelSize,vizStruct); 
camlight headlight; 

%%

imOrigin=min(V_grid,[],1)-voxelSize/2;

%%

controlPar_isosurface.nSub=[1 1 1];%round(max(v)/2./v);
controlPar_isosurface.capOpt=0; %Option to cap open ended surfaces
controlPar_isosurface.voxelSize=voxelSize;
controlPar_isosurface.contourLevel=contourLevel;
[Fi,Vi]=levelset2isosurface(M,controlPar_isosurface); %Get iso-surface
% Fi=fliplr(Fi); %Invert face orientation
[Fi,Vi]=mergeVertices(Fi,Vi); %merge nodes (iso-surface can be sloppy)
Fi_sorted=sort(Fi,2);
logicInvalid=any(diff(Fi_sorted,1,2)==0,2);
Fi=Fi(~logicInvalid,:);
[Fi,Vi]=patchCleanUnused(Fi,Vi);

Vi=Vi(:,[2 1 3]);
Vi=Vi+imOrigin;

% [Fi,Vi]=subtri(Fi,Vi,1);

[Fi,Vi]=triSurfRemoveThreeConnect(Fi,Vi);
[Fi,Vi]=patchCleanUnused(Fi,Vi);

%Smoothen angle of boundaries while holding radii constant

Eb=patchBoundary(Fi,Vi);
G=tesgroup(Eb); %group boundary sets

for q=1:1:size(G,2)
    logicGroupNow=G(:,q);
    
    indBoundaryNow=unique(Eb(logicGroupNow,:));
    logicBoundary=false(size(Vi,1),1);
    logicBoundary(indBoundaryNow)=1;

    Vi_centre=mean(Vi(indBoundaryNow,:),1);
    Vi=Vi-Vi_centre;
        
    [~,~,r] = cart2sph(Vi(:,1),Vi(:,2),Vi(:,3));        
    
    controlPar_smooth.Method='HC';
    controlPar_smooth.Alpha=0.1;
    controlPar_smooth.Beta=0.5;
    controlPar_smooth.n=25;
    controlPar_smooth.RigidConstraints=find(~logicBoundary);
    [Vi(:,[2 3])]=patchSmooth(Fi,Vi(:,[2 3]),[],controlPar_smooth);
    [az,elev,~] = cart2sph(Vi(:,1),Vi(:,2),Vi(:,3));
   
    [Vi(:,1),Vi(:,2),Vi(:,3)] = sph2cart(az,elev,r);
    Vi=Vi+Vi_centre;
end

%Smoothen surface while holding on to boundary   
Eb=patchBoundary(Fi,Vi);
controlPar_smooth.Method='HC';
controlPar_smooth.Alpha=0.1;
controlPar_smooth.Beta=0.5;
controlPar_smooth.n=150;
controlPar_smooth.RigidConstraints=unique(Eb(:));
[Vi]=patchSmooth(Fi,Vi,[],controlPar_smooth);

optionStruct.maxAngleDeviation=10;
optionStruct.triangleConvert=1;
[Fi,Vi]=tri2quadGroupSplit(Fi,Vi,optionStruct);

Eb=patchBoundary(Fi,Vi);
controlPar_smooth.Method='HC';
controlPar_smooth.Alpha=0.1;
controlPar_smooth.Beta=0.5;
controlPar_smooth.n=150;
controlPar_smooth.RigidConstraints=unique(Eb(:));
[Vi]=patchSmooth(Fi,Vi,[],controlPar_smooth);

%%

Vi=Vi-imOrigin; 
Vi=Vi(:,[2 1 3]);
Fi=fliplr(Fi); 
Eb=patchBoundary(Fi,Vi);

gpatch(Fi,Vi,'w','k',1,2);

%%

cFigure; 
gpatch(Fi,Vi,'w','k');
patchNormPlot(Fi,Vi);
axisGeom;
view(-90,0);
camlight headlight;
drawnow; 

%%

%Control settings
cPar.sphereRadius=sphereRadius;
cPar.coreRadius=cPar.sphereRadius/2;
cPar.numElementsMantel=numElementsMantel; 
cPar.numElementsCore=round(numElementsMantel*1.5); 
cPar.makeHollow=0;
cPar.cParSmooth.n=25;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
Es=meshStruct.E; %The elements 
Vs=meshStruct.V; %The vertices
Fb=meshStruct.Fb; %The boundary faces

%%
%Create cut view
Y=Vs(:,2); YE=mean(Y(Es),2);
L=YE>mean(Y);
[Fs,~]=element2patch(Es(L,:),[],'hex8');

cFigure;
subplot(1,2,1); hold on;
title('The hexahedral mesh sphere','FontSize',fontSize);
gpatch(Fb,Vs,'r');
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on;
title('Cut-view of the mesh','FontSize',fontSize);
gpatch(Fs,Vs,'r');
axisGeom(gca,fontSize);
camlight headlight;

drawnow; 

%% Shift branch set

Vi=Vi(:,[2 1 3]);
Vi=Vi+imOrigin;

% Eb=patchBoundary(Fi,Vi);
% G=tesgroup(Eb); %group boundary sets
% [~,indGet]=max(sum(G,1));
% logicGet=G(:,indGet);
% 
% [indGet]=edgeListToCurve(Eb(logicGet,:));
% indGet=indGet(1:end-1);
% 
% xc=mean(Vi(indGet,1));
% [yc,zc]=polycentroid(Vi(indGet,2)',Vi(indGet,3)');
%  
% Vc=[xc yc zc];
% % Vi=Vi-Vc; 
% 
% % Vi=Vi-mean(Vi(indGet,:));
% 
% cFigure; hold on; 
% gpatch(Fi,Vi,'w','k',1,3);
% plotV(Vi(indGet,:),'b-','LineWidth',3);
% axisGeom;
% camlight headlight;
% drawnow; 


%% Join model node sets

V=[Vs; Vi; ];
F_tube=Fi+size(Vs,1);
F_tube=fliplr(F_tube);

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb,V,'rw','k',0.8);
hl(2)=gpatch(F_tube,V,'bw','k',0.8);
patchNormPlot(F_tube,V);
legend(hl,{'Clot','Vessel'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Get contact surfaces
%

F_contact_blob=Fb;

%%
% Visualize contact surfaces

cFigure; hold on;
title('Tube blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','k',1);
patchNormPlot(F_tube,V);
hl(2)=gpatch(F_contact_blob,V,'gw','k',1);
patchNormPlot(F_contact_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis.ATTR.type=analysisType;
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=step_size;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;
febio_spec.Control.symmetric_stiffness=symmetric_stiffness; 
febio_spec.Control.min_residual=min_residual;

%Material section
if viscoElastic==1
    %Viscoelastic part
    febio_spec.Material.material{1}.ATTR.type='viscoelastic';
    febio_spec.Material.material{1}.ATTR.Name='Block_material';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.g1=g1;
    febio_spec.Material.material{1}.t1=t1;
    febio_spec.Material.material{1}.density=d;
    
    %Elastic part
    febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden unconstrained';
    febio_spec.Material.material{1}.elastic{1}.c1=c1;
    febio_spec.Material.material{1}.elastic{1}.m1=m1;
    febio_spec.Material.material{1}.elastic{1}.c2=c1;
    febio_spec.Material.material{1}.elastic{1}.m2=-m1;
    febio_spec.Material.material{1}.elastic{1}.cp=k;
    febio_spec.Material.material{1}.elastic{1}.density=d;
else
    febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.c1=c1;
    febio_spec.Material.material{1}.m1=m1;
    febio_spec.Material.material{1}.c2=c1;
    febio_spec.Material.material{1}.m2=-m1;
    febio_spec.Material.material{1}.cp=k;
    febio_spec.Material.material{1}.density=d;
end

febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=mean(Vi,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Clot'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(Es,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=Es;

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Vessel'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(Es,1)+(1:1:size(F_tube,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=F_tube;

% % -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(F_tube,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=F_tube;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave1';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_blob,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_blob;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%Boundary condition section 
febio_spec.Loads.body_load{1}.ATTR.type='const';
febio_spec.Loads.body_load{1}.x.VAL=bodyLoadMagnitude;
febio_spec.Loads.body_load{1}.x.ATTR.lc=1;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

%Contact section
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.two_pass=1;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1;
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

%LoadData
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; timeTotal 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='internal';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if 1%runFlag==1 %i.e. a succesful run
    
       
    %%
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    
    %%  
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    gpatch(F_tube,V,'w','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]/3);
    axis(axisLim(V_DEF)); %Set axis limits statically  
    camlight headlight;
    drawnow; 
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitue
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;

end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
