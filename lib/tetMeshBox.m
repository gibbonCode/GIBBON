function [meshStruct]=tetMeshBox(boxDim,pointSpacing)


%% Create surface mesh
[F,V,faceBoundaryMarker]=triBox(boxDim,pointSpacing);

%%  Mesh model using tetrahedral elements using tetGen

[regionA]=tetVolMeanEst(F,V); %Volume for regular tets

stringOpt='-pq1.2AaYQ';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=mean(V,1); %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker


[meshStruct]=runTetGen(inputStruct); %Run tetGen 


