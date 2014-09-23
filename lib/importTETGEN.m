function [meshOutput]=importTETGEN(loadNameStruct)

%%

dispStartTitleGibbonCode('Importing TetGen files');

%% IMPORT ELE file
fid=fopen(loadNameStruct.loadName_ele,'r');
[A]=textscan(fid,'%d %d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
elementID=double(A{1});
E=nan(max(elementID),4);
E(elementID,:)=double([A{2} A{3} A{4} A{5}]);
elementMaterialID=double(A{6});

%% IMPORT NODE file
fid=fopen(loadNameStruct.loadName_node,'r');
[A]=textscan(fid,'%d %f %f %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
nodeID=double(A{1});
V=nan(max(nodeID),3);
V(nodeID,:)=[A{2} A{3} A{4}];

%% IMPORT FACES file
fid=fopen(loadNameStruct.loadName_face,'r');
[A]=textscan(fid,'%d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
faceID=double(A{1});
F=nan(max(faceID),3);
F(faceID,:)=double([A{2} A{3} A{4}]);
faceBoundaryID=double(A{5});

%% CONVERT ELEMENTS TO FACES
[FE,faceMaterialID]=element2patch(E,elementMaterialID);

%% Create meshOutput structure
meshOutput.nodes=V; 
meshOutput.facesBoundary=F; 
meshOutput.boundaryMarker=faceBoundaryID;
meshOutput.faces=FE; 
meshOutput.elements=E; 
meshOutput.elementMaterialID=elementMaterialID; 
meshOutput.faceMaterialID=faceMaterialID; 

%%
disp(['--- Done --- ',datestr(now)]);

