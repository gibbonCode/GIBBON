function [meshOutput]=runTetGenSmesh(smeshStruct)

% function [meshOutput]=runTetGenSmesh(smeshStruct)
% ------------------------------------------------------------------------
% This function creates a .smesh files which it passes to TETGEN
% (<http://wias-berlin.de/software/tetgen/>) for tetrahedral meshing. The
% input structure should contain the following fields: 
% 
% smeshStruct.stringOpt=stringOpt;
% Where stringOpt should contain a strin composed of valid TETGEN command
% line switches (e.g. '-pq1.2VAaY'). See info below and the TETGEN manual.
%
% smeshStruct.Faces=F; 
% Where F is an array for all the model (triangular) faces
%
% smeshStruct.Nodes=V; 
% Where V is the array containing the nodal coordinates
%
% smeshStruct.holePoints=V_holes;
% Where V_holes describes a point (not part of V) that lies inside a
% prescribed hole inside the mesh. If V_holes=[] no holes are defined. 
%
% smeshStruct.faceBoundaryMarker=faceBoundaryMarker; 
% Where faceBoundaryMarker defines a label for each face denoting its
% membership to a particular boundary (e.g. all outer faces could have the
% same boundary label while an internal set of faces defining a hole has a
% different label). 
%
% smeshStruct.regionPoints=V_regions;
% Where similarly to V_holes the array V_regions defines points (not part
% of V) that lie inside a specific region (e.g. a material can be contained
% within another material and this allows you to specify multiple materials
% with different mesh densities and output labels).
%
% smeshStruct.regionA=regionA;
% Where regionA is a vector defining the A specification (volume) for the
% corrensponding regions (as defined in V_regions)
%
% smeshStruct.minRegionMarker=2; %Minimum region marker
% Arbitrary region marker label. Regions are labeled minRegionMarker:1:...
% for all regions. 
%
% smeshStruct.smeshName=smeshName; 
% Where smeshName is the file name for the input .smesh file. Only .smesh
% files are currently supported as input files. This function generates the
% .smesh file using the |writeBasicSmesh| function
%
%
% Below is a list of command line switches from the user manual:
% 
% -p Tetrahedralizes a piecewise linear complex (PLC).
% -Y Preserves the input surface mesh (does not modify it).
% -r Reconstructs a previously generated mesh.
% -q Refines mesh (to improve mesh quality).
% -R Mesh coarsening (to reduce the mesh elements).
% -A Assigns attributes to tetrahedra in different regions.
% -a Applies a maximum tetrahedron volume constraint.
% -m Applies a mesh sizing function.
% -i Inserts a list of additional points.
% -O Specifies the level of mesh optimization.
% -S Specifies maximum number of added points.
% -T Sets a tolerance for coplanar test (default 10?8).
% -X Suppresses use of exact arithmetic.
% -M No merge of coplanar facets or very close vertices.
% -w Generates weighted Delaunay (regular) triangulation.
% -c Retains the convex hull of the PLC.
% -d Detects self-intersections of facets of the PLC.
% -z Numbers all output items starting from zero.
% -f Outputs all faces to .face file.
% -e Outputs all edges to .edge file.
% -n Outputs tetrahedra neighbors to .neigh file.
% -v Outputs Voronoi diagram to files.
% -g Outputs mesh to .mesh file for viewing by Medit.
% -k Outputs mesh to .vtk file for viewing by Paraview.
% -J No jettison of unused vertices from output .node file.
% -B Suppresses output of boundary information.
% -N Suppresses output of .node file.
% -E Suppresses output of .ele file.
% -F Suppresses output of .face and .edge file.
% -I Suppresses mesh iteration numbers.
% -C Checks the consistency of the final mesh.
% -Q Quiet: No terminal output except errors.
% -V Verbose: Detailed information, more terminal output.
% -h Help: A brief instruction for using TetGen.
%
% See the TETGEN manual for more information.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/03/05
%------------------------------------------------------------------------

%%

dispStartTitleGibbonCode('TETGEN Tetrahedral meshing');

%% PARSE INPUT STRUCTURE
stringOpt=smeshStruct.stringOpt;
smeshName=smeshStruct.smeshName;

%% PARSE FILE NAMES
[modelPath,modelFileName,~]=fileparts(smeshName);
modelName=fullfile(modelPath,modelFileName);
loadName_ele=[modelName,'.1.ele'];
loadName_node=[modelName,'.1.node'];
loadName_face=[modelName,'.1.face'];
loadName_edge=[modelName,'.1.edge'];

loadNameStruct.loadName_ele=loadName_ele;
loadNameStruct.loadName_node=loadName_node;
loadNameStruct.loadName_face=loadName_face;
loadNameStruct.loadName_edge=loadName_edge;

%% BUILD SMESH FILE
writeBasicSmesh(smeshStruct);

%% SETTING TETGEN PATHNAMES
% try
%     pathNameTetGen='C:\Users\kmmoerman\00_WORK\SOURCE_CODES\tetGen\tetgen1.5.0\Release';
%     pathNameTetView='C:\Users\kmmoerman\00_WORK\SOURCE_CODES\tetGen\tetView';
% catch
    pathNameTetGen=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen');
    pathNameTetView=pathNameTetGen;
% end
runNameTetGen=fullfile(pathNameTetGen,'tetgen.exe');

%% DELETE POSSIBLE EXISTING TETGEN FILES

if exist(loadNameStruct.loadName_ele,'file')
    delete(loadNameStruct.loadName_ele);
end

if exist(loadNameStruct.loadName_node,'file')    
    delete(loadNameStruct.loadName_node);
end

if exist(loadNameStruct.loadName_face,'file')
    delete(loadNameStruct.loadName_face);
end

%% RUN TETGEN
disp(['--- Running TetGen for meshing --- ',datestr(now)]);
runString=[runNameTetGen,' ',stringOpt,' ',smeshName];
system(runString);
dispDoneGibbonCode;

% %% COPY OUTPUT FILES TO TETVIEW DIRECTORY
% smeshName_copy=fullfile(pathNameTetView,[modelFileName,'.smesh']);
% copyfile(smeshName,smeshName_copy);
% loadName_ele_copy=fullfile(pathNameTetView,[modelFileName,'.1.ele']);
% copyfile(loadName_ele,loadName_ele_copy);
% loadName_node_copy=fullfile(pathNameTetView,[modelFileName,'.1.node']);
% copyfile(loadName_node,loadName_node_copy);
% loadName_face_copy=fullfile(pathNameTetView,[modelFileName,'.1.face']);
% copyfile(loadName_face,loadName_face_copy);
% loadName_edge_copy=fullfile(pathNameTetView,[modelFileName,'.1.edge']);
% copyfile(loadName_edge,loadName_edge_copy);

%% IMPORT TETGEN FILES

try
    [meshOutput]=importTETGEN(loadNameStruct);
    meshOutput.loadNameStruct=loadNameStruct;
catch
   error('TetGen output not found!');
end

