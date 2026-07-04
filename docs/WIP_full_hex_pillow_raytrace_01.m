%% WIP_full_hex_pillow_raytrace_01
% This code may evolve into functionality or a demo to show full hexahedral
% meshing using "pillowing" and ray-tracing. 

%% Desription
% This code aims to create hexahedral meshes on arbetrary surfaces. It
% first converts the surface to a 3D image. Next voxels fully inside the
% surface are identified. Next this "interior image" is eroded such that is
% a certain distance away from the input surface geometry. The eroded voxel
% set is used as a hexahedral mesh. Next this hexahedral mesh is "pillowed"
% outward towards the input surface. This outward offsetting is based on
% ray-tracing of the interior hexahedral mesh boundary surface towards the
% encapsulating input mesh. During ray-tracing the normal directions for
% ray-tracing, as well as the ray-traced output mesh, are smoothed and the
% ray tracing can be repeated iteratively to let the mesh "settle" on the
% input mesh surface. 
%
% TO DO: 
% * Check robustness of ray-tracing implementation
% * Improve efficiency, perhaps use C++ ray-tracing code
% * Check for optimal conditions for the erosion/growth of the voxel mesh
% in terms of pillowed mesh quality. 
% * Check for optimal smoothing settings

%% 

clear; close all; clc;

%% Plot settings 
fontSize=15;
faceAlpha1=1;
faceAlpha2=0.3;
edgeColor=0.3*ones(1,3);
edgeWidth=1; 
plotColors=gjet(4);

%% Control parameters

%Ray tracing settings
rayTraceOptionStruct.tolEps         = 1e-6;
rayTraceOptionStruct.triSide        = -1;
rayTraceOptionStruct.rayType        = 'line';
rayTraceOptionStruct.exclusionType  = 'inclusive';
rayTraceOptionStruct.paired         = 0;

nErode=2; %Number of erosion steps

numSmoothQuad=500; %Number of smoothing steps of interior boundary surface
numSmoothQuadNormals=25; %Number of normal smoothing steps prior to ray tracing 
numSmoothHex=25; %Number of smoothing steps for hexahedral mesh (while boundary is held fixed). 

numStepsRayUpdate=6; %Number of ray-tracing repetitions
numPillowSteps=nErode; %Number of elements in pillowed mesh

%% Create test geometry 

testCase=4; 
switch testCase
    case 1
        [Fs,Vs]=geoSphere(3,1);
        voxelSize=0.1;
    case 2
        [Fs,Vs]=graphicsModels(5);
        Vs=Vs*1000;
        [Fs,Vs]=triSurfRemoveThreeConnect(Fs,Vs);
        voxelSize=2; % The output image voxel size.
    case 3
        [Fs,Vs]=stanford_bunny;
        voxelSize=2;
    case 4
        %Torus parameters
        r=1; %Sphere radius
        rc=2.5; %Central radius
        n=3;
        nr=16*n;
        nc=25*n;
        patchType='tri';
        
        [Fs,Vs]=patchTorus(r,nr,rc,nc,patchType);
        voxelSize=0.2;
end

%% Create surface image 

% Using |triSurf2Im| function to convert patch data to image data
[M,G,~]=triSurf2Im(Fs,Vs,voxelSize);
imOrigin=G.origin;
L_model=(M==2); %Interior only choosen here

%% 
% Visualize surface image

sliceViewerOption.edgeColor='k';
sv3(M,voxelSize,sliceViewerOption);
gpatch(Fs,Vs-imOrigin,'g','none',faceAlpha2);
[~,hc]=icolorbar;
set(hc,'YTickLabel',{'Exterior','Boundary','Intertior'});
drawnow;

%% Erode to form interior mesh image
 
%Defining erosion/dilation kernel
adjacencyMask=sphereIm(1,1);
for q=1:1:nErode
    L_model = convn(double(L_model),adjacencyMask,'same')==7;  %Erode
end

%Remove poorly connected voxels
L_model_ini=convn(double(L_model),adjacencyMask,'same')>=2;
while 1
    nn=nnz(L_model);
    L_model = convn(double(L_model),adjacencyMask,'same')>=3 & L_model_ini;
    if nnz(L_model)==nn
        break
    end
end

%% Visualize interior mesh image

sliceViewerOption.edgeColor='k';
sv3(L_model,voxelSize,sliceViewerOption);
gpatch(Fs,Vs-imOrigin,'g','none',faceAlpha2);
% icolorbar;
drawnow;

%% Get hexahedral elements for interior voxels
% 

[E_hex,V_hex,C_hex]=im2patch(M,L_model,'h');
[F_hex,C_hex_F]=element2patch(E_hex,C_hex); % Element faces
% [F_hex,~,~]=uniqueIntegerRow(F_hex); %Unique face set (only 1 of each shared face is kept)
[F_hex,V_hex,~,indFix]=mergeVertices(F_hex,V_hex);
E_hex=indFix(E_hex);

% Convert Coordinates
[V_hex(:,1),V_hex(:,2),V_hex(:,3)]=im2cart(V_hex(:,2),V_hex(:,1),V_hex(:,3),voxelSize*ones(1,3));
V_hex=V_hex+imOrigin(ones(size(V_hex,1),1),:);

%Get mesh boundary faces
indBoundaryFaces=tesBoundary(F_hex);
Fb=F_hex(indBoundaryFaces,:);
indBoundaryNodes=unique(Fb(:));

%Isolate boundary face mesh from rest
[Fq,Vq,indFix]=patchCleanUnused(Fb,V_hex);
indBoundary_q=indFix(indBoundaryNodes);

%% Visualize interior hex mesh

cFigure;
title('Visualizing pillowing process','FontSize',fontSize);
hold on;

gpatch(Fb,V_hex,'gw','k',1,1);
hp=gpatch(Fq,Vq,'rw','k',0.75,1);
gpatch(Fs,Vs,'w','none',0.25);
% patchNormPlot(Fq,Vq)
% patchNormPlot(Fs,Vs)
axisGeom(gca,fontSize); camlight headlight
gdrawnow;

%% Use ray-tracing and mesh pillowing to form outer mesh layer

cParNorm.Method='HC';
cParNorm.n=numSmoothQuadNormals;

cParSurf.Method='HC';
cParSurf.n=numSmoothQuad;

Vqr=Vq;
for q=1:1:numStepsRayUpdate   
    if q==1
        dMax=voxelSize*10;
    else
        dMax=voxelSize*2;
    end

    [Vqr]=patchSmooth(Fq,Vqr,[],cParSurf);    
    [~,~,Nq]=patchNormal(Fq,Vqr);
    [Nq]=patchSmooth(Fq,Nq,[],cParNorm);
    
    [Vqr]=rayTraceOut(Vqr-((voxelSize/100)*Nq),Nq,Fs,Vs,rayTraceOptionStruct,dMax);
       
    set(hp,'Vertices',Vqr);
    drawnow;
end

%% Compose hex layer of pillowed elements

[Ep,Vp]=fromToMesh(Fq,Vq,Vqr,numPillowSteps) ;
[Fp,Cp]=element2patch(Ep,[]); % Element faces

%% Join sets and merge nodes

[E,V,C]=joinElementSets({E_hex,Ep},{V_hex,Vp});
F=element2patch(E,C,'hex8');
[F,V,~,indFix]=mergeVertices(F,V);
E=indFix(E);

%% Retrieve boundary faces and nodes
[indBoundaryFaces]=freeBoundaryPatch(F);
Fb=F(indBoundaryFaces,:);
indBoundaryNodes=unique(Fb(:));

%% Smoothen hex mesh 
% During smoothing boundary mesh is held fixed

cParHex.Method='HC';
cParHex.n=numSmoothHex;
cParHex.RigidConstraints=indBoundaryNodes;
[V]=tesSmooth(E,V,[],cParHex);

%% Visualizing final mesh

meshStruct.elements=E;
meshStruct.nodes=V;
meshStruct.faces=F;
meshStruct.elementMaterialID=C;
meshViewOptionStruct.cutDir=2;
meshViewOptionStruct.numSliceSteps=50;
meshViewOptionStruct.edgeColor='w';
meshView(meshStruct,meshViewOptionStruct);

%%


[A,EE,AE]=dihedralAngles(E,V,'hex8');
A=180*(A./pi);
AE=180*(AE./pi);

A_max=max(A,[],2);
A_min=min(A,[],2);

[F,A_max_F]=element2patch(E,A_max);
[~,A_min_F]=element2patch(E,A_min);

%%

cFigure; 
subplot(1,2,1); hold on;
title(['Max dihedral angle ',num2str(max(A_max_F))])
gpatch(F,V,A_max_F,'k',1,1);
axisGeom; camlight headlight; 
colormap(gca,gjet(25)); colorbar; 
clim([min(A(:)) max(A(:))]);

subplot(1,2,2); hold on;
title(['Min dihedral angle ',num2str(min(A_min_F))])
gpatch(F,V,A_min_F,'k',1,1);
axisGeom; camlight headlight; 
colormap(gca,gjet(25)); colorbar; 
clim([min(A(:)) max(A(:))]);
gdrawnow; 


%%

function [V1r]=rayTraceOut(V1,N1,F2,V2,optStruct,dMax)

V1r=V1;
numSteps=size(V1,1);
hw=waitbar(0,['Ray tracing...',100,'%']);
for q=1:1:numSteps
    v1=V1(q,:);
    n1=N1(q,:);
    [V_intersect,~,d_intersect]=triSurfRayTrace(v1,dMax*n1,F2,V2,optStruct);
    if ~isempty(V_intersect)
        d=abs(d_intersect);
        [~,indMin]=min(d);
        V_now=V_intersect(indMin,:);
        V1r(q,:)=V_now;
    end

    waitbar(q/numSteps,hw,['Ray tracing...',num2str(round(100.*q/numSteps)),'%']);
end
close(hw);
end

%%

function [varargout]=fromToMesh(Fp1,Vp1,Vp2,numSteps) 

%Get coordinates
X=linspacen(Vp1(:,1),Vp2(:,1),numSteps+1);
Y=linspacen(Vp1(:,2),Vp2(:,2),numSteps+1);
Z=linspacen(Vp1(:,3),Vp2(:,3),numSteps+1);

%Collect node set
V=[X(:) Y(:) Z(:)];

%% Create element sets 

if isa(Fp1,'cell')
    E=Fp1;
    Fp1n=Fp1;
    Fp2n=Fp1;
    for q=1:1:numel(Fp1)
        [E{q},Fp1n{q},Fp2n{q}]=createElementSet(Fp1{q},Vp1,numSteps);
    end
else     
    [E,Fp1n,Fp2n]=createElementSet(Fp1,Vp1,numSteps);    
end

%% Collect output
varargout{1}=E;
varargout{2}=V;
varargout{3}=Fp1n;
varargout{4}=Fp2n;
 
end

function [E,Fp1,Fp2]=createElementSet(Fp1,Vp1,numSteps)

%Create element matrix
E=repmat(Fp1,[numSteps,2]);
E_add=0:size(Vp1,1):size(Vp1,1)*(numSteps-1);
E_add=E_add(ones(size(Fp1,1),1),:);
E_add=E_add(:);
E_add=E_add(:,ones(size(Fp1,2),1));
E_add=[E_add E_add+size(Vp1,1)];
E=E+E_add;

%Create top and bottom face set
Fp1=E(1:size(Fp1,1),1:size(Fp1,2));
Fp2=E(1+(end-size(Fp1,1)):end,size(Fp1,2)+1:end);

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
