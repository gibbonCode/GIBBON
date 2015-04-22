function [meshOutput]=importTETGEN(loadNameStruct)

%%

dispStartTitleGibbonCode('Importing TetGen files');

%% IMPORT ELE file
try
    [~,E,elementMaterialID]=importEleFile_tetGen(loadNameStruct.loadName_ele);
catch
%     disp([loadNameStruct.loadName_ele,' import unsuccesful']);
    E=[];
    elementMaterialID=[];
end

%% IMPORT NODE file
try
    [~,V]=importNodeFile_tetGen(loadNameStruct.loadName_node);
catch
    warning([loadNameStruct.loadName_node,' import unsuccesful']);
    V=[];
end

%% IMPORT FACES file
try
    [~,F,faceBoundaryID]=importFaceFile_tetGen(loadNameStruct.loadName_face);
catch
%     warning([loadNameStruct.loadName_face,' import unsuccesful']);
    F=[];
    faceBoundaryID=[];
end

%% CONVERT ELEMENTS TO FACES
if ~isempty(E)
    [FE,faceMaterialID]=element2patch(E,elementMaterialID,'tet4');
else
    FE=[];
    faceMaterialID=[];
end

%% Create meshOutput structure
meshOutput.nodes=V; 
meshOutput.facesBoundary=F; 
meshOutput.boundaryMarker=faceBoundaryID;
meshOutput.faces=FE; 
meshOutput.elements=E; 
meshOutput.elementMaterialID=elementMaterialID; 
meshOutput.faceMaterialID=faceMaterialID; 
meshOutput.loadNameStruct=loadNameStruct; 

%%
disp(['--- Done --- ',datestr(now)]);

