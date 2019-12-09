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
F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};
[F,V]=mergeVertices(F,V);

R=euler2DCM([0 0 0.5*pi]);
V=V*R;

%%
cFigure; 
gpatch(F,V,'kw')
axisGeom;
camlight headlight;
drawnow;

%%
[regionA]=tetVolMeanEst(F,V); %Volume for regular tets
stringOpt='-pq1.2AaYQ';
inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=ones(size(F,1),1); %Face boundary markers
inputStruct.regionPoints=[0 0 0]; %region points
inputStruct.regionA=regionA*5;

[meshStruct]=runTetGen(inputStruct); %Run tetGen

%%


E=meshStruct.elements;
Fb=meshStruct.facesBoundary;
V=meshStruct.nodes;

%%
% Visualize input geometry 

meshView(meshStruct);

%% Example: Create a dual lattice mesh with outer surface cladding
% See also: |dualClad|

t=linspace(min(V(:,2))+0.5,max(V(:,2))-0.5,12);
for q=1:1:12
    logicLattice_V=V(:,2)<t(q);
    logicLattice_E=all(logicLattice_V(E),2);
    
    FT_other_all=element2patch(E(~logicLattice_E,:),V);
    indB=tesBoundary(FT_other_all,V);
    FT_other=FT_other_all(indB,:);
    
    cladOpt=1;
    shrinkFactor=0.2;
    [FT,VT,CT]=dualLattice(E(logicLattice_E,:),V,shrinkFactor,cladOpt);
    
    %%
    % Visualize results
    
    cFigure; hold on;
    gpatch(FT_other,V,'w','none',1);
    gpatch(FT,VT,CT,'none',1);
    % patchNormPlot(FT,VT);
    colormap(viridis(4)); icolorbar;
    axis equal; axis vis3d;
    view(2);
    camlight headlight;
    axis off; colorbar off;
    drawnow;
    
    fileName=['/mnt/data/MATLAB/temp/efw/hipImplantMeshLattice_',num2str(q),'.png'];
    export_fig(fileName,'-png','-r100');
    
end
