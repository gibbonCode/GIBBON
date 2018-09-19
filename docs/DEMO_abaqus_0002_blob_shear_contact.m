%% DEMO_febio_0031_blob_shear_contact
% Below is a demonstration for:
% 
% * Building geometry for a hemi-spherical blob with tetrahedral elements
% which is being sheared by a rigid wall. 
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
% * contact, sliding, sticky, friction
% * rigid body constraints
% * tetrahedral elements, tet4
% * triangular elements, tri3
% * slab, block, rectangular
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

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
abaqusInpFileNamePart='tempModel';
abaqusInpFileName=fullfile(savePath,[abaqusInpFileNamePart,'.inp']); %INP file name

% Hemi-sphere parameters
hemiSphereRadius=1; 
nRefine=1; 
closeOption=1; 
smoothEdge=1; 

% Ground plate parameters
plateRadius=2*hemiSphereRadius; 

% Probe parameters
probeWidth=3*hemiSphereRadius; 
filletProbe=0.5; %Fillet radius

% Define probe displacement
probeDisplacement=hemiSphereRadius*2; 
proveOverlapFactor=0.4;

% Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=15;
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=25; %Maximum number of retires
symmetric_stiffness=0;
min_residual=1e-20;
step_size=1/numTimeSteps;
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/(numTimeSteps); %Maximum time step size

%Contact parameters
contactPenalty=20;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.1; 

%% Creating model geometry and mesh
% 

% Create hemi-sphere mesh
[F_blob,V_blob,C_blob]=hemiSphereMesh(nRefine,hemiSphereRadius,1);
pointSpacingBlob=mean(patchEdgeLengths(F_blob,V_blob));

%Smoothen edges
if smoothEdge==1
    %Get rigid region
    ind=1:1:size(V_blob,1); %Indices for all nodes
    indRigid=find(ismember(ind,F_blob(C_blob==2,:)) & ~ismember(ind,F_blob(C_blob==1,:))); %Indices for new bottom surface nodes
    
    %Smoothing
    cPar.Method='HC';
    cPar.n=250;
    cPar.RigidConstraints=indRigid;
    [V_blob]=patchSmooth(F_blob,V_blob,[],cPar);
    
    %Fix color data with new bottom surface
    C_blob=ones(size(C_blob));
    C_blob(all(ismember(F_blob,indRigid),2))=2;

end

%%
% Visualize hemi-sphere surface

cFigure; hold on;
gtitle('The hemi-sphere surface mesh',fontSize);
gpatch(F_blob,V_blob,C_blob,'k',0.85);
patchNormPlot(F_blob,V_blob);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 
drawnow; 

%% 
% Using tetgen to create a tetrahedral mesh

% Tetgen input structure
inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=F_blob;
inputStruct.Nodes=V_blob;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=C_blob; %Face boundary markers
inputStruct.regionPoints=getInnerPoint(F_blob,V_blob); %region points
inputStruct.regionA=tetVolMeanEst(F_blob,V_blob);
inputStruct.minRegionMarker=2; %Minimum region marker
 
% Create tetrahedral mesh using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 
 
% Access model element and patch data
Fb_blob=fliplr(meshOutput.facesBoundary);
Cb_blob=meshOutput.boundaryMarker;
V_blob=meshOutput.nodes;
E_blob=meshOutput.elements;

%%
% Visualize blob mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
gpatch(Fb_blob,V_blob,Cb_blob,'k',1);
patchNormPlot(Fb_blob,V_blob);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb_blob,V_blob,'kw','none',0.25);
meshView(meshOutput,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

%% Creating rigid body ground plate
% 

%Get outer surve of ground surface 
[Eb]=patchBoundary(Fb_blob(Cb_blob==2,:),V_blob);
indCurveBottom=edgeListToCurve(Eb);
indCurveBottom=indCurveBottom(1:end-1);

% Derive point spacing for plate
pointSpacingPlate=pointSpacingBlob; 

% Compose outer curve of the plate
nPlateCurve=ceil((2*pi*plateRadius)/pointSpacingPlate);
t=linspace(0,2*pi,nPlateCurve);
t=t(1:end-1); 
x=plateRadius.*sin(t);
y=plateRadius.*cos(t); 
Vp_outer_curve=[x(:) y(:)];

% Copy inner curve from the hemi-sphere
Vp_inner_curve=V_blob(indCurveBottom,[1 2]);

% Create mesh out outer region
regionCell={Vp_outer_curve,Vp_inner_curve};
[F_plate,V_plate]=regionTriMesh2D(regionCell,pointSpacingPlate,0,0);
V_plate(:,3)=0; %Add z-direction

% Copy mesh for inner region
[F_plate_inner,V_plate_inner]=patchCleanUnused(fliplr(Fb_blob(Cb_blob==2,:)),V_blob);
[F_plate,V_plate,C_plate]=joinElementSets({F_plate,F_plate_inner},{V_plate,V_plate_inner});
[F_plate,V_plate]=mergeVertices(F_plate,V_plate);

center_of_mass_plate=mean(V_plate,1);

%%
% Visualizing plate mesh

cFigure; hold on;
gtitle('The plate surface mesh',fontSize);
gpatch(Fb_blob,V_blob,'kw','none',0.5);
gpatch(F_plate,V_plate,C_plate,'k',1);
patchNormPlot(F_plate,V_plate);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 
drawnow; 

%% Creating rigid body shear surface

pointSpacingProbe=pointSpacingBlob/2; 

%Sketching side profile
x=[-hemiSphereRadius hemiSphereRadius hemiSphereRadius]-hemiSphereRadius*2;
y=[0 0 0];
z=[hemiSphereRadius*(1-proveOverlapFactor) hemiSphereRadius*(1-proveOverlapFactor) hemiSphereRadius*1.5];
V_probe_curve_sketch=[x(:) y(:) z(:)];

%Fillet sketch
np=100; %Number of points used to construct each fillet edge
[V_probe_curve]=filletCurve(V_probe_curve_sketch,filletProbe,np,0);
numPointsProbeCurve=ceil(max(pathLength(V_probe_curve))/pointSpacingProbe);
[V_probe_curve] = evenlySampleCurve(V_probe_curve,numPointsProbeCurve,'pchip',0);

% Extruding curve
% controlParametersExtrude.pointSpacing=pointSpacingProbe;
controlParametersExtrude.depth=hemiSphereRadius*2.5; 
controlParametersExtrude.numSteps=ceil(controlParametersExtrude.depth/pointSpacingProbe);
controlParametersExtrude.numSteps=controlParametersExtrude.numSteps+iseven(controlParametersExtrude.numSteps); %Force uneven
controlParametersExtrude.patchType='tri'; 
controlParametersExtrude.dir=0;
controlParametersExtrude.n=[0 1 0];
controlParametersExtrude.closeLoopOpt=0; 

[F_probe,V_probe]=polyExtrude(V_probe_curve,controlParametersExtrude);
F_probe=fliplr(F_probe); %Invert face orientation so normals point to blob

center_of_mass_probe=mean(V_probe,1);

%%
% Visualizing probe mesh

cFigure; hold on;
title('The probe surface mesh','fontSize',fontSize);
gpatch(Fb_blob,V_blob,'kw','none',0.5);
gpatch(F_plate,V_plate,'kw','none',0.5);
hl(1)=plotV(V_probe_curve_sketch,'k.-.','lineWidth',3,'MarkerSize',25);
hl(2)=plotV(V_probe_curve,'r-','lineWidth',3,'MarkerSize',25);
hl(3)=gpatch(F_probe,V_probe,'gw','k',1);
legend(hl,{'Sketched probe curve','Rounded probe curve','Probe surface mesh'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Join model node sets

V=[V_blob; V_plate; V_probe];
F_plate=F_plate+size(V_blob,1);
F_probe=F_probe+size(V_blob,1)+size(V_plate,1);
Fb_all=[Fb_blob;F_plate;F_probe];

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
hl(2)=gpatch(F_plate,V,'bw','k',0.8);
hl(3)=gpatch(F_probe,V,'gw','k',0.8);
legend(hl,{'Blob','Plate','Probe'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Get contact surfaces
%

% F_contact_blob1=Fb_blob(Cb_blob==1,:);
% F_contact_blob2=Fb_blob(Cb_blob==2,:);

%%
% Visualize contact surfaces

cFigure; 
subplot(1,2,1); hold on;
title('Probe blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_probe,V,'rw','k',1);
patchNormPlot(F_probe,V);
hl(2)=gpatch(Fb_blob,V,'gw','k',1);
patchNormPlot(Fb_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 

subplot(1,2,2); hold on;
title('Plate blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_plate,V,'rw','k',1);
patchNormPlot(F_plate,V);
hl(2)=gpatch(Fb_blob,V,'gw','k',1);
patchNormPlot(Fb_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 

drawnow; 

%% Defining the abaqus input structure
% See also |abaqusStructTemplate| and |abaqusStruct2inp| and the abaqus user
% manual.

%%--> Heading
abaqus_spec.Heading.COMMENT{1}='Job name: ABAQUS inp file creation demo';
abaqus_spec.Heading.COMMENT{2}='Generated by: GIBBON';

%%--> Preprint
abaqus_spec.Preprint.ATTR.echo='NO';
abaqus_spec.Preprint.ATTR.model='NO';
abaqus_spec.Preprint.ATTR.history='NO';
abaqus_spec.Preprint.ATTR.contact='NO';

%--> Part

% Node
nodeIds=(1:1:size(V,1))';
abaqus_spec.Part.COMMENT='This section defines the part geometry in terms of nodes and elements';
abaqus_spec.Part.ATTR.name='Cube';
abaqus_spec.Part.Node={nodeIds,V};

% Element
elementIds=(1:1:size(E_blob,1))';
abaqus_spec.Part.Element{1}.ATTR.type='C3D8';%'C3D8R';
abaqus_spec.Part.Element{1}.VAL={elementIds,E_blob};

% Element sets
abaqus_spec.Part.Elset{1}.ATTR.elset='Set-1';
abaqus_spec.Part.Elset{1}.VAL=elementIds;

% Sections
abaqus_spec.Part.Solid_section.ATTR.elset='Set-1';
abaqus_spec.Part.Solid_section.ATTR.material='Elastic';

%%--> Assembly
abaqus_spec.Assembly.ATTR.name='Assembly-1';
abaqus_spec.Assembly.Instance.ATTR.name='Test-assembly';
abaqus_spec.Assembly.Instance.ATTR.part='Test';

abaqus_spec.Assembly.Nset{1}.ATTR.nset='All';
abaqus_spec.Assembly.Nset{1}.ATTR.instance=abaqus_spec.Assembly.Instance.ATTR.name;
abaqus_spec.Assembly.Nset{1}.VAL=[1:1:size(V,1)];

%%--> Material
abaqus_spec.Material.ATTR.name='Elastic';
abaqus_spec.Material.Elastic=[0.5 0.4];

%%--> Step
abaqus_spec.Step.ATTR.name='Step-1';
abaqus_spec.Step.ATTR.nlgeom='YES';
abaqus_spec.Step.Static=[0.1 1 1e-5 0.1];

% Boundary
% abaqus_spec.Step.Boundary{1}.VAL={'Set-1',[1,1]};
% abaqus_spec.Step.Boundary{2}.VAL={'Set-2',[2,2]};
% abaqus_spec.Step.Boundary{3}.VAL={'Set-3',[3,3]};
% abaqus_spec.Step.Boundary{4}.VAL={'Set-4',[3,3],-0.1};

%Output
abaqus_spec.Step.Restart.ATTR.write='';
abaqus_spec.Step.Restart.ATTR.frequency=0;

abaqus_spec.Step.Output{1}.ATTR.field='';
abaqus_spec.Step.Output{1}.ATTR.variable='PRESELECT';
abaqus_spec.Step.Output{2}.ATTR.history='';
abaqus_spec.Step.Output{2}.ATTR.variable='PRESELECT';
% abaqus_spec.Step.Node_print.ATTR.nset='all';
% abaqus_spec.Step.Node_print.ATTR.frequency = 1;
% abaqus_spec.Step.Node_print.VAL='COORD';
% abaqus_spec.Step.El_print.VAL='S';

%%

[T]=abaqusStruct2inp(abaqus_spec,abaqusInpFileName);

% textView(abaqusInpFileName);

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