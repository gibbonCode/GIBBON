function [T]=writeBasicSmesh(smeshStruct)

%%

dispStartTitleGibbonCode('Writing SMESH file');

%% PARSE INPUT STRUCTURE

F=smeshStruct.Faces;
V=smeshStruct.Nodes;
V_holes=smeshStruct.holePoints;
faceBoundaryMarker=smeshStruct.faceBoundaryMarker;
V_regions=smeshStruct.regionPoints; %region points
regionA=smeshStruct.regionA;
minRegionMarker=smeshStruct.minRegionMarker; %Minimum region marker
smeshName=smeshStruct.smeshName; 

%% PART 1 NODES
disp('----> Adding node field');

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

%% PART 2 FACETS
disp('----> Adding facet field');

F_field=[size(F,2).*ones(size(F,1),1) F faceBoundaryMarker(:)];

%Create text form for facet field
textForm=repmat('%i ',1,2+size(F,2)); %Integer parts
textForm=textForm(1:end-1); %Remove last space
textForm(end+1:end+2)='\n'; %Add end of line part
F_char=sprintf(textForm,F_field');
F_cell=regexp(F_char,'\n','split')';

if numel(F_cell)>1
    F_cell=F_cell(1:end-1);
end

T(end+1)={'#PART 2 - Facet list'};
T(end+1)={'#num faces, boundary markers'};
numFacets=size(F,1);
boundMarker=1;
doubleList=[numFacets boundMarker];
charList=sprintf('%i %i',doubleList');
T(end+1)={charList};
T(end+1)={'#Facet ID, <corner1, corner2, corner3,...>,[attribute],[boundary marker]'};
T(end+1:end+numel(F_cell))=F_cell;

%% PART 3 HOLES
disp('----> Adding holes specification');

if ~isempty(V_holes) %If holes are present
    V_id=1:size(V_holes,1);
    V_field=[V_id(:) V_holes];
    V_char=sprintf('%i %0.16e %0.16e %0.16e \n',V_field');
    V_cell=regexp(V_char,' \n','split')';
    
    T(end+1)={'#PART 3 - Hole list'};
    T(end+1)={'#Num holes'};
    T(end+1)={sprintf('%i',size(V_holes,1))};
    T(end+1)={'#<hole #> <x> <y> <z>'};
    T(end+1:end+numel(V_cell))=V_cell;
else %No holes present
    T(end+1)={'#PART 3 - Hole list'};
    T(end+1)={'0'};
end

%% PART 4 REGIONS
disp('----> Adding region specification');

regionMarkers=minRegionMarker:1:size(V_regions,1)+minRegionMarker-1;

V_field=[regionMarkers(:) V_regions -regionMarkers(:) regionA(:)];
V_char=sprintf('%i %0.16e %0.16e %0.16e %i %0.16e \n',V_field');
V_cell=regexp(V_char,' \n','split')';

T(end+1)={'#PART 4 - Region list'};
T(end+1)={'#Num regions'};
T(end+1)={sprintf('%i',numel(regionMarkers))};
T(end+1)={'#<region #> <x> <y> <z> <region number> <region attribute>'};
T(end+1:end+numel(V_cell))=V_cell;

%% SAVING TXT FILE
cell2txtfile(smeshName,T,0);

dispDoneGibbonCode;

