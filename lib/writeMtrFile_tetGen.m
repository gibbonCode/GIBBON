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



 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
