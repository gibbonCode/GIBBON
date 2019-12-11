clear; close all; clc; 

%%

markerSize1=25; 
fontSize=15; 

%%

% Load surface geometry
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 

stlName='hip_implant_iso_merge.stl';
fileName=fullfile(pathName,stlName); 
[stlStruct] = import_STL(fileName);
F1=stlStruct.solidFaces{1};
V1=stlStruct.solidVertices{1};
[F1,V1]=mergeVertices(F1,V1);

Q=eye(3,3);
V1=V1*10;
R=euler2DCM([0 -0.5*pi 0]);
V1=V1*R;
Q=Q*R;
R=euler2DCM([0 0 pi]);
V1=V1*R;
Q=Q*R;
R=euler2DCM([-0.215*pi 0 0]);
V1=V1*R;
Q=Q*R;
R=euler2DCM([0 -0.05*pi 0]);
V1=V1*R;
Q=Q*R;
R=euler2DCM([0 0 0.15*pi]);
V1=V1*R;
Q=Q*R;
V1=V1*Q';

stlName='femur_iso.stl';
fileName=fullfile(pathName,stlName); 
[stlStruct] = import_STL(fileName);
F2=stlStruct.solidFaces{1};
V2=stlStruct.solidVertices{1};
[F2,V2]=mergeVertices(F2,V2);
V2=V2*1000;
R=euler2DCM([-0.5*pi 0 0.5*pi]);
V2=V2*R;
R=euler2DCM([-0.05*pi 0 0]);
V2=V2*R;
R=euler2DCM([0 0 0.33*pi]);
V2=V2*R;
V2=V2*Q';

%%

cFigure; hold on;
gpatch(F1,V1,'kw','k',1);
gpatch(F2,V2,'bw','k',0.5);
axisGeom;
camlight headlight;
drawnow; 

%%
D=minDist(V2,V1);
logicCut=V2(:,1)<=49 & V2(:,2)>=min(V1(:,2))-5 & D<25;

logicCut=all(logicCut(F2),2);
logicCut=triSurfLogicSharpFix(F2,logicCut,3);

cFigure; hold on;
gpatch(F1,V1,'kw','k',1);
gpatch(F2,V2,logicCut,'k',0.5);
axisGeom;
camlight headlight;
drawnow; 

%%
[F2,V2]=patchCleanUnused(F2(~logicCut,:),V2);

Eb2=patchBoundary(F2,V2);
ind2=edgeListToCurve(Eb2);
ind2=ind2(1:end-1);

%%

logicCut1=V1(:,1)<=43;
logicCut1=all(logicCut1(F1),2);
logicCut1=triSurfLogicSharpFix(F1,logicCut1,3);

Eb1=patchBoundary(F1(logicCut1,:),V1);
ind1=edgeListToCurve(Eb1);
ind1=ind1(1:end-1);

[f,v]=regionTriMesh3D({V2(ind2,:),V1(ind1,:)},[],'linear');
[F2,V2,C2]=joinElementSets({F2,f},{V2,v});
[F2,V2]=mergeVertices(F2,V2);

%%
cFigure; hold on;
gpatch(F1,V1,'kw','k',1);
gpatch(F2,V2,C2,'k',1);

axisGeom;
camlight headlight;
drawnow; 

%%

logic2=C2==2; 
for q=1:1:2
    ind2=unique(F2(logic2,:));
    logic2=any(ismember(F2,ind2),2);
end
indRigid=unique(F2(~logic2,:));

Eb=patchBoundary(F2,V2);
indRigid=unique([indRigid;Eb(:)]);

cPar.n=5;
cPar.Method='HC';
cPar.RigidConstraints=indRigid;
[V2]=patchSmooth(F2,V2,[],cPar);

cFigure; hold on;
gpatch(F1,V1,logicCut1,'k',1);
gpatch(F2,V2,'kw','none',0.5);

axisGeom;
camlight headlight;
drawnow; 

%%


Eb=patchBoundary(F1(logicCut1,:),V1);
ind=unique(Eb(:));
[D2,indMin]=minDist(V2,V1(ind,:));

Eb2=patchBoundary(F2,V2);
indMove=find(D2<10);
indMove=indMove(~ismember(indMove,Eb2));

[N,Vn,Nv]=patchNormal(F1,V1);

d=mean(patchEdgeLengths(F2,V2));

V2(indMove,:)=V2(indMove,:)+0.25*d*Nv(ind(indMin(indMove)),:);

%%
cFigure; hold on;

gpatch(F1,V1,'kw','k',1);
patchNormPlot(F1,V1);
gpatch(F2,V2,D2,'none',0.5);
plotV(V2(indMove,:),'r.','markerSize',markerSize1);

axisGeom;
camlight headlight;
drawnow; 

% 
% ind1=unique(F1(logicCut1,:));
% D=minDist(V2,V1(ind1,:))<5;
% 
% cFigure; hold on;
% gpatch(F1,V1,'kw','none',0.5);
% gpatch(F2,V2,D,'k',1);
% 
% axisGeom;
% camlight headlight;
% drawnow; 

%%

[FT,VT,CT]=joinElementSets({F1,F2},{V1,V2},{double(logicCut1),2*ones(size(F2,1),1)});
[FT,VT]=mergeVertices(FT,VT);
VT=VT*Q;

cFigure; hold on;
gpatch(FT,VT,CT,'k',0.5);

axisGeom;
camlight headlight;
colormap gjet; icolorbar;
drawnow; 

%%



%%

%%
% Define region points
[V_region1]=getInnerPoint(FT(CT==0 | CT==2,:),VT);
[V_region2]=getInnerPoint(FT(CT==0 | CT==1,:),VT);

%%
% Visualize interior points

cFigure; hold on;
gpatch(FT,VT,'kw','none',0.2);
plotV(V_region1,'r.','markerSize',markerSize1);
plotV(V_region2,'b.','markerSize',markerSize1);
camlight('headlight'); 
axisGeom(gca,fontSize);
drawnow; 

%% 
% Mesh using tetgen
inputStruct.stringOpt='-pq1.2AaY'; %TetGen option string
inputStruct.Faces=FT; %The faces
inputStruct.Nodes=VT; %The vertices
inputStruct.holePoints=[]; %The hole interior points
inputStruct.faceBoundaryMarker=CT; %Face boundary markers
inputStruct.regionPoints=[V_region1;V_region2]; %The region interior points
inputStruct.regionA=[tetVolMeanEst(F2,V2) tetVolMeanEst(F1,V1)]*4; %Volume for regular tets

%% 
% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb_foot=meshOutput.facesBoundary; %Boundary faces of the foot
Cb_foot=meshOutput.boundaryMarker; %Boundary marker/color data for the foot
V_foot=meshOutput.nodes; %The vertices/nodes
E_foot=meshOutput.elements; %The tet4 elements

%% 
% Visualizing mesh using |meshView|, see also |anim8|
optionStruct.cutDir=1;
meshView(meshOutput,optionStruct);

%%
% [F1,V1]=patchCleanUnused(F(G==1,:),V);
% [F2,V2]=patchCleanUnused(F(G==2,:),V);
% 
% N2=patchNormal(F2,V2);
% n=[-1 0 0];
% D=dot(N2,n(ones(size(N2,1),1),:),2);
% logicLeft=D>0.9;
% logicLeft=triSurfLogicSharpFix(F2,logicLeft,3);
% 
% [F2,V2]=patchCleanUnused(F2(~logicLeft,:),V2);
% 
% N1=patchNormal(F1,V1);
% n=[1 0 0];
% D=dot(N1,n(ones(size(N1,1),1),:),2);
% logicLeft=D>0.999;
% logicLeft=triSurfLogicSharpFix(F1,logicLeft,3);
% [F1,V1]=patchCleanUnused(F1(~logicLeft,:),V1);
% 
% Eb1=patchBoundary(F1,V1);
% ind1=edgeListToCurve(Eb1);
% ind1=ind1(1:end-1);
% Eb2=patchBoundary(F2,V2);
% ind2=edgeListToCurve(Eb2);
% ind2=ind2(1:end-1);
% ind2=flip(ind2);
% 
% %Create input structure
% [f,v,boundaryInd]=regionTriMesh3D({V1(ind1,:),V2(ind2,:)},[],'linear');
% 
% %%
% 
% cFigure; hold on;
% gpatch(F1,V1,'rw','none',0.25);
% gpatch(F2,V2,'bw');
% gpatch(f,v,'gw');
% plotV(V1(ind1,:),'r-','LineWidth',3);
% plotV(V2(ind2,:),'b-','LineWidth',3);
% plotV(V1(ind1(1),:),'r.','MarkerSize',30);
% plotV(V2(ind2(1),:),'b.','MarkerSize',30);
% axisGeom;
% camlight headlight;
% colorbar;
% drawnow; 
% 
% %%
% 
% [F,V,C]=joinElementSets({F1,F2,f},{V1,V2,v});
% 
% [F,V]=mergeVertices(F,V);
% 
% cFigure; hold on;
% gpatch(F,V,C,'k',1);
% 
% axisGeom;
% camlight headlight;
% colorbar;
% drawnow; 
% 
% stlStruct.solidNames={'hip_model'};
% stlStruct.solidVertices={V};
% stlStruct.solidFaces={F};
% fileName=fullfile(pathName,'hip_implant_iso_fix.stl'); 
% export_STL_txt(fileName,stlStruct);
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
