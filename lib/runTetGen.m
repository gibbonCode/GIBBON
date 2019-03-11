function [meshOutput]=runTetGen(inputStruct)

% function [meshOutput]=runTetGen(smeshStruct)
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
% Change log: 
% 2014/03/05 Created as runTetGenSmesh
% 2014/10/20 Copied and renamed as runTetGen
% 2014/10/20 Altered to allow for constrained and unconstrained Delaunay
% tesselation (rather than meshing). This means the input file can be a
% .node file instead of a smesh file. These changes are made to allow for
% sizing function specification on the initial delaunay mesh.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/20
% 2015/06/23 Fixed tet10 and boundaryMarker handling
% 2015/07/10 Fixed handling of spaces in paths
% 2015/09/22 Updated for new tet4_tet10 function, removed for loop for
% boundaryMarker handling for tet10 elements
% 2018/05/15 Change temp files to be the general GIBBON temp folder
%------------------------------------------------------------------------

%% PARSE INPUT

if ~isfield(inputStruct,'tetType')
    inputStruct.tetType='tet4'; %Revert to default if missing
elseif isempty(inputStruct.tetType)
    inputStruct.tetType='tet4'; %Revert to default if missing
end

switch inputStruct.tetType
    case 'tet4' %Linear tetrahedral elements

    case 'tet10' %Quadratic tetrahedral elements

    otherwise
        error('Wrong element type specified, valid options are tet4 or tet10, for other element types use converions (e.g. tet2hex)');
end

%%

dispStartTitleGibbonCode('TETGEN Tetrahedral meshing');

%% PARSE INPUT STRUCTURE

if ~isfield(inputStruct,'stringOpt')
    inputStruct.stringOpt='-pY';
elseif isempty(inputStruct.stringOpt)
    inputStruct.stringOpt='';
end

if ~isfield(inputStruct,'Faces')
    inputStruct.Faces=[];
end

if isempty(inputStruct.Faces)
    useDelaunayConstraints=0;
else
    useDelaunayConstraints=1;
end

if ~isfield(inputStruct,'Nodes')
    error('No nodes specified in input structure');
end

if ~isfield(inputStruct,'holePoints')
    inputStruct.holePoints=[];
end

if ~isfield(inputStruct,'faceBoundaryMarker')
    inputStruct.faceBoundaryMarker=ones(size(inputStruct.Faces,1),1);
end

if ~isfield(inputStruct,'regionPoints')
    inputStruct.regionPoints=[];
end

if ~isfield(inputStruct,'regionA')
    inputStruct.regionA=[];
end

if ~isfield(inputStruct,'minRegionMarker')
    inputStruct.minRegionMarker=2;
end

if ~isfield(inputStruct,'modelName')
    if isfield(inputStruct,'smeshName') %WILL BE REMOVED
        inputStruct.modelName=inputStruct.smeshName;
        warning('smeshStruct.smeshName input will be replaced by smeshStruct.modelName in future releases!');    
        copyFiles=1;
    else
        inputStruct.modelName=[]; 
        copyFiles=0;
    end
else
    copyFiles=1;
end

if isfield(inputStruct,'sizeData')    
    sizingOn=1; 
else
    sizingOn=0; 
end

%% PARSE FILE NAMES

%Setting tetGen pathnames
pathNameTetGen=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen');
pathNameTempFiles=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','temp');

%Create temp folder if it does not exist
if ~exist(pathNameTempFiles,'file')
    mkdir(pathNameTempFiles);
end

%Get model name
modelName=inputStruct.modelName;
if isempty(modelName)
    modelName=fullfile(pathNameTempFiles,'temp');
end

%Remove extension if present
[savePathStr,modelNameClean,~] = fileparts(modelName);
modelName=fullfile(savePathStr,modelNameClean);
modelNameTemp=fullfile(pathNameTempFiles,modelNameClean);

inputStruct.modelName=modelNameTemp;
if isfield(inputStruct,'smeshName') %WILL BE REMOVED
    inputStruct.smeshName=modelNameTemp;
end
% inputStruct = rmfield(inputStruct,'smeshName');

cleanUpTetGen(pathNameTempFiles); % Clean up temp directory

%Create possible output file names
loadNameStructTemp.loadName_ele =[modelNameTemp,'.1.ele'];
loadNameStructTemp.loadName_node=[modelNameTemp,'.1.node'];
loadNameStructTemp.loadName_face=[modelNameTemp,'.1.face'];
loadNameStructTemp.loadName_edge=[modelNameTemp,'.1.edge'];
loadNameStructTemp.loadName_smesh=[modelNameTemp,'.smesh'];

loadNameStruct.loadName_ele =[modelName,'.1.ele'];
loadNameStruct.loadName_node=[modelName,'.1.node'];
loadNameStruct.loadName_face=[modelName,'.1.face'];
loadNameStruct.loadName_edge=[modelName,'.1.edge'];
loadNameStruct.loadName_smesh=[modelName,'.smesh'];

%% Set tetgen run patch string depending on operational system

compString=computer; 
switch compString
    case 'PCWIN' %Windows 32-bit
        error('PCWIN 32-bit is not supported. Compile tetGen from the source and alter the code here');
    case 'PCWIN64' %Windows 64-bit
        pathNameTetGenFile=fullfile(pathNameTetGen,'win64');
        runNameTetGen=fullfile(pathNameTetGenFile,'tetgen.exe');
    case 'GLNXA64'        
        pathNameTetGenFile=fullfile(pathNameTetGen,'lin64');
        runNameTetGen=fullfile(pathNameTetGenFile,'tetgen');
    case 'MACI64'        
        pathNameTetGenFile=fullfile(pathNameTetGen,'mac64');
        runNameTetGen=fullfile(pathNameTetGenFile,'tetgen');
    otherwise
        error('Your platform does not seem to be supported');
end

%% Create input file which is either a .node or a .smesh file

switch useDelaunayConstraints
    case 0
        writeNodeFile_tetGen(inputStruct);
        runModelName=[modelNameTemp,'.node'];
    case 1
        writeBasicSmesh(inputStruct);
        runModelName=[modelNameTemp,'.smesh'];
end

%% RUN TETGEN
bOpt=0;
if sizingOn    
    
    indNodes=unique(inputStruct.Faces(:));
    numNodes=numel(indNodes);
    
    if numNodes<size(inputStruct.Nodes,1) %If the number of nodes is larger than those used in the boundary
        bOpt=1;
        
        %Compute Delauney tesselation first
        
        %Pass on only Q option if provided
        if strfind(inputStruct.stringOpt,'Q')
            strOpt='-Q';
        else
            strOpt='';
        end
        
        runString=['"',runNameTetGen,'" ',strOpt,' "',runModelName,'"']; % Old which may not work with paths with spaces:  runString=[runNameTetGen,' ',strOpt,' ',runModelName];
        
        disp(['--- Running TetGen to compute Delaunay tesselation of input set --- ',datestr(now)]);
        [runStatus,runOut]=system(runString,'-echo');
        dispDoneGibbonCode;
        
        %Rename .1.node files to .b.node
        copyfile(loadNameStructTemp.loadName_node,[modelNameTemp,'.b.node']);
        
        %Rename .1.ele files to .b.ele
        copyfile(loadNameStructTemp.loadName_ele,[modelNameTemp,'.b.ele']);
        
        %Delete existing files
        if exist(loadNameStructTemp.loadName_ele,'file')==2
            delete(loadNameStructTemp.loadName_ele);
        end
        
        if exist(loadNameStructTemp.loadName_node,'file')==2
            delete(loadNameStructTemp.loadName_node);
        end
        
        if exist(loadNameStructTemp.loadName_face,'file')==2
            delete(loadNameStructTemp.loadName_face);
        end
        
    end

    %Delete existing .mtr files
    extCell={'mtr'}; %Extensions of files to delete
    
    for qc=1:1:numel(extCell)
        ext=extCell{qc}; %Current extension
        fileList = dir(fullfile(pathNameTempFiles,['*.',ext]));
        fileList={fileList(1:end).name}; %Current file list
        
        %Delete files
        for q=1:1:numel(fileList)
            fileName=fullfile(savePathStr,fileList{q});
            delete(fileName);
        end
    end
        
    %Create .mtr file
    writeMtrFile_tetGen(inputStruct,bOpt);
    
    %Compute Delauney tesselation / Mesh using tetgen

    runString=['"',runNameTetGen,'" ',[inputStruct.stringOpt,'m'],' "',runModelName,'"']; 
            
    disp(['--- Running TetGen to mesh initial Delaunay tesselation using sizing function--- ',datestr(now)]);
    [runStatus,runOut]=system(runString,'-echo');
    dispDoneGibbonCode;
else
    %Compute Delauney tesselation / Mesh using tetgen
    runString=['"',runNameTetGen,'" ',inputStruct.stringOpt,' "',runModelName,'"']; 
    disp(['--- Running TetGen to mesh input boundary--- ',datestr(now)]);
    [runStatus,runOut]=system(runString,'-echo');
    dispDoneGibbonCode;
end

%% IMPORT TETGEN FILES

try
    [meshOutput]=importTETGEN(loadNameStructTemp); 
    meshOutput.loadNameStruct=loadNameStruct; 
catch
   error('Error importing TetGen output files. Check file names.');
end


%% Convert element type if required

switch inputStruct.tetType
    case 'tet4' %Linear tetrahedral elements
        %Keep as is
    case 'tet10' %Quadratic tetrahedral elements
        
        E_tet4=meshOutput.elements;
        V_tet4=meshOutput.nodes;        
        
        Fb_tet4=meshOutput.facesBoundary;
               
        
        elementMaterialID=meshOutput.elementMaterialID;    
        
        [E_tet10,V_tet10,~,Fb_tet10,~]=tet4_tet10(E_tet4,V_tet4,[],Fb_tet4);

        [F_tet10,faceMaterialID]=element2patch(E_tet10,elementMaterialID,'tet10');   

        % Compose output          
        meshOutput.nodes=V_tet10;
        meshOutput.facesBoundary=Fb_tet10;
%         meshOutput.boundaryMarker=faceBoundaryID; %Remains valid
        meshOutput.faces=F_tet10;
        meshOutput.elements=E_tet10;        
%         meshOutput.elementMaterialID=elementMaterialID; %Remains valid
        meshOutput.faceMaterialID=faceMaterialID;     
          
end

%% Copy relevant files

if copyFiles
    extCell={'ele','node','face','edge','smesh'}; %Extensions of files to copy    
    for qc=1:1:numel(extCell)
        ext=extCell{qc}; %Current extension
        fileList = dir(fullfile(pathNameTempFiles,['*.',ext]));
        fileList={fileList(1:end).name}; %Current file list
        
        %Copying files to output location
        for q=1:1:numel(fileList)
            fileNameTemp=fullfile(pathNameTempFiles,fileList{q});            
            fileName=fullfile(savePathStr,fileList{q});            
            if ~strcmp(fileNameTemp,fileName) %If the names are not the same
                %A location other than the temp folder is an output folder
                copyfile(fileNameTemp,fileName);
            end
        end
    end    
end

%% Clean up directory
cleanUpTetGen(pathNameTempFiles); % Clean up temp directory
 
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
