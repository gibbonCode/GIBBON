function writeMtrFile_tetGen(inputStruct,bOpt)

%%

disp(['--- Writing MTR file --- ',datestr(now)]); %Start message

%% PARSE INPUT STRUCTURE

sizeData=inputStruct.sizeData;
mtrFileName=inputStruct.modelName; 

%Force extension to be .node
[pathstr,name,~] = fileparts(mtrFileName);
if bOpt
    mtrFileName=fullfile(pathstr,[name,'.b.mtr']);
else
    mtrFileName=fullfile(pathstr,[name,'.mtr']);
end

%%


V_field=sizeData(:);
V_char=sprintf('%0.16e \n',V_field');
V_cell = regexp(V_char, '\n', 'split')'; 

if numel(V_cell)>1
    V_cell=V_cell(1:end-1);
end

T=cell(1,1);
T(1,1)={'#num nodes, attribute size'};
numNodes=numel(sizeData);
attributeSize=1; 

doubleList=[numNodes attributeSize];
charList=sprintf('%i %i ',doubleList');
T(end+1,1)={charList};
T(end+1,1)={'#attribute'};
T(end+1:end+numel(V_cell),1)=V_cell;

%% SAVING TXT FILE
cell2txtfile(mtrFileName,T,0);



