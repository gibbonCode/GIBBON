%% WIP_DEMO_gmsh_meshing
% Below is a demonstration for:
% 
% * Building geometry for a slab and dome using gmsh
% * Importing and visualizing the displacement results

%% Keywords
%
% * gmsh
%%
clear; close all; clc;


%%
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
%path with pre-written .geo (gmsh) files
geo_path = fullfile(defaultFolder,'data','gmsh');
%list of prewritten .geo files used to generate meshes
geo_files = {'\rigid_dome.geo','\square_membrane.geo'};
%path to gmsh generated mesh files
mesh_path = savePath;
%list of gmsh generated mesh files
mesh_files = {'\rigid_dome.msh','\square_membrane.msh'};
%path to an installation of gmsh
gmsh_path = fullfile(defaultFolder,'lib_ext','gmsh','win64');

%%

parameter_names = {'membrane_L', 'membrane_W', 'membrane_H', 'radius'};
parameter_values = {50, 50, 0.5, 12.5};

%update the geo files 
for geo_ind = 1:length(geo_files)

    geo_text = txtfile2cell([geo_path,geo_files{geo_ind}]);
    
    for param_ind = 1:length(parameter_names)
        parameter_line =  find(startsWith(geo_text,parameter_names{param_ind}));
        if ~isempty(parameter_line)
            geo_text{parameter_line} = [parameter_names{param_ind}, ' = ', num2str(parameter_values{param_ind}),';'];
            
        end
    end
    
    cell2txtfile([savePath,geo_files{geo_ind}],geo_text);

end

%% Load Meshes
for i = 1:size(geo_files,2)
    runGmsh([savePath,geo_files{i}],gmsh_path)
end


%% Import Meshes
for i = 1:size(mesh_files,2)

Meshes{i} = read_gmsh([mesh_path, mesh_files{i}]);

end

%% Combine Meshes
MeshStruct = Meshes{1};
if size(Meshes,2)>1
    for i = 2:size(Meshes,2)

    vol_ct = max(MeshStruct.elementMaterialID);
    surf_ct = max(MeshStruct.boundaryMarker);
    
    MeshStruct_temp = Meshes{i};
    
    nodes = MeshStruct_temp.nodes;
    node_IDs = MeshStruct_temp.node_IDs+max(MeshStruct.node_IDs);
    facesBoundary = MeshStruct_temp.facesBoundary+max(MeshStruct.node_IDs);
    boundaryMarker = MeshStruct_temp.boundaryMarker+surf_ct;
    elements = MeshStruct_temp.elements+max(MeshStruct.node_IDs);
    elementMaterialID = MeshStruct_temp.elementMaterialID+vol_ct;
    part_name = MeshStruct_temp.loadNameStruct.MeshName;
    volumes_names = MeshStruct_temp.loadNameStruct.VolumeNames;
    facesBoundary_Names = MeshStruct_temp.loadNameStruct.SurfaceNames;
    

    MeshStruct.nodes = [MeshStruct.nodes;nodes];
    MeshStruct.facesBoundary = [MeshStruct.facesBoundary;facesBoundary];
    MeshStruct.boundaryMarker = [MeshStruct.boundaryMarker;boundaryMarker];
    MeshStruct.elements = [MeshStruct.elements;elements];
    MeshStruct.elementMaterialID = [MeshStruct.elementMaterialID;elementMaterialID];
    MeshStruct.loadNameStruct.MeshName = {MeshStruct.loadNameStruct.MeshName;part_name};
    MeshStruct.loadNameStruct.VolumeNames = [MeshStruct.loadNameStruct.VolumeNames;volumes_names];
    MeshStruct.loadNameStruct.SurfaceNames = [MeshStruct.loadNameStruct.SurfaceNames;facesBoundary_Names];
    

    end
end

vol_ct = max(MeshStruct.elementMaterialID);
surf_ct = max(MeshStruct.boundaryMarker);


%% 
% % Plotting model boundary surfaces and a cut view
% 
E1=MeshStruct.elements; %The elements 
V1=MeshStruct.nodes; %The nodes (vertices)
Fb1=MeshStruct.facesBoundary; %The boundary faces
Cb1=MeshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
%elementMaterialIndices=ones(size(E1,1),1); %Element material indices

%% Display the mesh

%Plot settings
fontSize=15;
faceAlpha1=1;
faceAlpha2=0.3;
markerSize=40;
markerSize2=20;
lineWidth=3;

hFig=cFigure;
hs1=subplot(1,2,1);
hold on;
for i = 1:surf_ct
    
    hl(i) = gpatch(Fb1(Cb1==i,:),V1,Cb1(Cb1==i),'k',faceAlpha1);

end
colormap(gjet(6))
legend(hl,MeshStruct.loadNameStruct.SurfaceNames,'Interpreter','none');
axisGeom(gca,fontSize);
drawnow;

hs2=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs2];
meshView(MeshStruct,optionStruct);
axisGeom(gca,fontSize);
C = findall(hFig,'type','ColorBar');
set(C,'TickLabels',MeshStruct.loadNameStruct.VolumeNames)
set(C,'TickLabelInterpreter','none')
%        'TickLabels',MeshStruct.loadNameStruct.VolumeNames,'TickLabelInterpreter','none')

drawnow;
%%



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
