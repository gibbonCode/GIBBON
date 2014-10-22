function [TR]=constrainedDelaunayTetGen(V,C)


%% Create TetGen input structure

if ~isempty(C)
    inputStruct.Faces=C; %Add face constraints if not empty
    inputStruct.stringOpt='-pQY';
else
    inputStruct.stringOpt='-Q'; 
end
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=ones(size(C,1),1); %Face boundary markers
inputStruct.regionPoints=[]; %region points

%% Run TetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

TR = triangulation(meshOutput.elements,meshOutput.nodes);

