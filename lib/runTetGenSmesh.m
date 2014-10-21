function [meshOutput]=runTetGenSmesh(inputStruct)

warning('runTetGenSmesh has been renamed to runTetGen. runTetGenSmesh will be removed in future releases.')
[meshOutput]=runTetGen(inputStruct); %Run tetGen

% 
% % function [meshOutput]=runTetGenSmesh(smeshStruct)
% % ------------------------------------------------------------------------
% % This function creates a .smesh files which it passes to TETGEN
% % (<http://wias-berlin.de/software/tetgen/>) for tetrahedral meshing. The
% % input structure should contain the following fields: 
% % 
% % smeshStruct.stringOpt=stringOpt;
% % Where stringOpt should contain a strin composed of valid TETGEN command
% % line switches (e.g. '-pq1.2VAaY'). See info below and the TETGEN manual.
% %
% % smeshStruct.Faces=F; 
% % Where F is an array for all the model (triangular) faces
% %
% % smeshStruct.Nodes=V; 
% % Where V is the array containing the nodal coordinates
% %
% % smeshStruct.holePoints=V_holes;
% % Where V_holes describes a point (not part of V) that lies inside a
% % prescribed hole inside the mesh. If V_holes=[] no holes are defined. 
% %
% % smeshStruct.faceBoundaryMarker=faceBoundaryMarker; 
% % Where faceBoundaryMarker defines a label for each face denoting its
% % membership to a particular boundary (e.g. all outer faces could have the
% % same boundary label while an internal set of faces defining a hole has a
% % different label). 
% %
% % smeshStruct.regionPoints=V_regions;
% % Where similarly to V_holes the array V_regions defines points (not part
% % of V) that lie inside a specific region (e.g. a material can be contained
% % within another material and this allows you to specify multiple materials
% % with different mesh densities and output labels).
% %
% % smeshStruct.regionA=regionA;
% % Where regionA is a vector defining the A specification (volume) for the
% % corrensponding regions (as defined in V_regions)
% %
% % smeshStruct.minRegionMarker=2; %Minimum region marker
% % Arbitrary region marker label. Regions are labeled minRegionMarker:1:...
% % for all regions. 
% %
% % smeshStruct.smeshName=smeshName; 
% % Where smeshName is the file name for the input .smesh file. Only .smesh
% % files are currently supported as input files. This function generates the
% % .smesh file using the |writeBasicSmesh| function
% %
% %
% % Below is a list of command line switches from the user manual:
% % 
% % -p Tetrahedralizes a piecewise linear complex (PLC).
% % -Y Preserves the input surface mesh (does not modify it).
% % -r Reconstructs a previously generated mesh.
% % -q Refines mesh (to improve mesh quality).
% % -R Mesh coarsening (to reduce the mesh elements).
% % -A Assigns attributes to tetrahedra in different regions.
% % -a Applies a maximum tetrahedron volume constraint.
% % -m Applies a mesh sizing function.
% % -i Inserts a list of additional points.
% % -O Specifies the level of mesh optimization.
% % -S Specifies maximum number of added points.
% % -T Sets a tolerance for coplanar test (default 10?8).
% % -X Suppresses use of exact arithmetic.
% % -M No merge of coplanar facets or very close vertices.
% % -w Generates weighted Delaunay (regular) triangulation.
% % -c Retains the convex hull of the PLC.
% % -d Detects self-intersections of facets of the PLC.
% % -z Numbers all output items starting from zero.
% % -f Outputs all faces to .face file.
% % -e Outputs all edges to .edge file.
% % -n Outputs tetrahedra neighbors to .neigh file.
% % -v Outputs Voronoi diagram to files.
% % -g Outputs mesh to .mesh file for viewing by Medit.
% % -k Outputs mesh to .vtk file for viewing by Paraview.
% % -J No jettison of unused vertices from output .node file.
% % -B Suppresses output of boundary information.
% % -N Suppresses output of .node file.
% % -E Suppresses output of .ele file.
% % -F Suppresses output of .face and .edge file.
% % -I Suppresses mesh iteration numbers.
% % -C Checks the consistency of the final mesh.
% % -Q Quiet: No terminal output except errors.
% % -V Verbose: Detailed information, more terminal output.
% % -h Help: A brief instruction for using TetGen.
% %
% % See the TETGEN manual for more information.
% %
% %
% % Kevin Mattheus Moerman
% % kevinmoerman@hotmail.com
% % 2014/03/05
% %------------------------------------------------------------------------
% 
% %% PARSE INPUT
% 
% if ~isfield(inputStruct,'tetType')
%     inputStruct.tetType='tet4'; %Revert to default if missing
% elseif isempty(inputStruct.tetType)
%     inputStruct.tetType='tet4'; %Revert to default if missing
% end
% 
% switch inputStruct.tetType
%     case 'tet4' %Linear tetrahedral elements
% 
%     case 'tet10' %Quadratic tetrahedral elements
% 
%     otherwise
%         error('Wrong element type specified, valid options are tet4 or tet10, for other element types use converions (e.g. tet2hex)');
% end
% 
% %%
% 
% dispStartTitleGibbonCode('TETGEN Tetrahedral meshing');
% 
% %% PARSE INPUT STRUCTURE
% stringOpt=inputStruct.stringOpt;
% 
% if isfield(inputStruct,'smeshName'); %WILL BE REMOVED
%     modelName=inputStruct.smeshName;
%     warning('smeshStruct.smeshName input will be replaced by smeshStruct.modelName in future releases!');
% elseif isfield(inputStruct,'modelName'); 
%     modelName=inputStruct.modelName; 
% end
% 
% %Remove possible extenion
% [pathstr,name,~] = fileparts(modelName);
% modelName=fullfile(pathstr,name);
% inputStruct.modelName=modelName;
% 
% %% PARSE FILE NAMES
% loadName_ele=[modelName,'.1.ele'];
% loadName_node=[modelName,'.1.node'];
% loadName_face=[modelName,'.1.face'];
% loadName_edge=[modelName,'.1.edge'];
% 
% loadNameStruct.loadName_ele=loadName_ele;
% loadNameStruct.loadName_node=loadName_node;
% loadNameStruct.loadName_face=loadName_face;
% loadNameStruct.loadName_edge=loadName_edge;
% 
% %% BUILD SMESH FILE
% 
% writeBasicSmesh(inputStruct);
% 
% %% SETTING TETGEN PATHNAMES
% 
% compString=computer; 
% switch compString
%     case 'PCWIN' %Windows 32-bit
%         error('PCWIN 32-bit is not supported. Compile tetGen from the source and alter the code here');
%     case 'PCWIN64' %Windows 64-bit
%         pathNameTetGen=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','win64');
%         runNameTetGen=fullfile(pathNameTetGen,'tetgen.exe');
%     case 'GLNXA64'        
%         pathNameTetGen=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','lin64');
%         runNameTetGen=fullfile(pathNameTetGen,'tetgen');
%     case 'MACI64'        
%         pathNameTetGen=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','mac64');
%         runNameTetGen=fullfile(pathNameTetGen,'tetgen');
%     otherwise
%         error('Your platform does not seem to be supported. Code your own solution or contact support.')
% end
% 
% %% DELETE POSSIBLE EXISTING TETGEN FILES
% 
% if exist(loadNameStruct.loadName_ele,'file')
%     delete(loadNameStruct.loadName_ele);
% end
% 
% if exist(loadNameStruct.loadName_node,'file')    
%     delete(loadNameStruct.loadName_node);
% end
% 
% if exist(loadNameStruct.loadName_face,'file')
%     delete(loadNameStruct.loadName_face);
% end
% 
% %% RUN TETGEN
% 
% smeshName=[modelName,'.smesh'];
% 
% disp(['--- Running TetGen for meshing --- ',datestr(now)]);
% runString=[runNameTetGen,' ',stringOpt,' ',smeshName];
% system(runString);
% dispDoneGibbonCode;
% 
% %% IMPORT TETGEN FILES
% 
% try
%     [meshOutput]=importTETGEN(loadNameStruct); 
% %     meshOutput.nodes=V;
% %     meshOutput.facesBoundary=F;
% %     meshOutput.boundaryMarker=faceBoundaryID;
% %     meshOutput.faces=FE;
% %     meshOutput.elements=E;
% %     meshOutput.elementMaterialID=elementMaterialID;
% %     meshOutput.faceMaterialID=faceMaterialID;
% %     meshOutput.loadNameStruct=loadNameStruct;
% catch
%    error('TetGen output not found!');
% end
% 
% %% Convert element type if required
% 
% switch inputStruct.tetType
%     case 'tet4' %Linear tetrahedral elements
%         %Keep as is
%     case 'tet10' %Quadratic tetrahedral elements
%         TET4=meshOutput.elements;
%         V4=meshOutput.nodes;
%         elementMaterialID=meshOutput.elementMaterialID;        
%         [TET10,V10,~]=tet4_tet10(TET4,V4,{});        
%         [F10,faceMaterialID]=element2patch(TET10,elementMaterialID,'tet10');   
%         matLabels=unique(elementMaterialID(:));
%         Fb=[];
%         faceBoundaryID=[];
%         for q=1:1:numel(matLabels)
%             matLabel=matLabels(q);
%             [Fbn,~,Cb]=meshBoundary(TET10(elementMaterialID==matLabel,:),'tet10',elementMaterialID(elementMaterialID==matLabel));
%             faceBoundaryIDn=Cb;
%             Fb=[Fb;Fbn];
%             faceBoundaryID=[faceBoundaryID;faceBoundaryIDn];            
%         end
%         
%         meshOutput.nodes=V10;
%         meshOutput.facesBoundary=Fb;
%         meshOutput.boundaryMarker=faceBoundaryID; 
%         meshOutput.faces=F10;
%         meshOutput.elements=TET10;
% %         meshOutput.elementMaterialID=elementMaterialID; %Remains valid
%         meshOutput.faceMaterialID=faceMaterialID;     
% end
% 
% 
