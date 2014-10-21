function [T]=writeNodeFile_tetGen(inputStruct)

%%

dispStartTitleGibbonCode('Writing NODE file');

%% PARSE INPUT STRUCTURE

V=inputStruct.Nodes;
nodeFileName=inputStruct.modelName; 

%Force extension to be .node
[pathstr,name,~] = fileparts(nodeFileName);
nodeFileName=fullfile(pathstr,[name,'.node']);

%%

V_id=1:1:size(V,1);
V_field=[V_id(:) V];
V_char=sprintf('%i %0.16e %0.16e %0.16e \n',V_field');
V_cell = regexp(V_char, '\n', 'split'); 

if numel(V_cell)>1
    V_cell=V_cell(1:end-1);
end

T={'#PART 1 - Node list';'#num nodes, num dimensions, num attributes, num boundary markers'};
numNodes=size(V,1);
numDims=3; 
numAtr=0;
boundMarker=0;
doubleList=[numNodes numDims numAtr boundMarker];
charList=sprintf('%i %i %i %i',doubleList');
T(end+1)={charList};
T(end+1)={'#Node ID, x, y, z,attribute,boundary marker'};
T(end+1:end+numel(V_cell))=V_cell;

%% SAVING TXT FILE

cell2txtfile(nodeFileName,T,0);
dispDoneGibbonCode;


