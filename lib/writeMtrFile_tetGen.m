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
cell2txtfile(mtrFileName,T,0,0);
 
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
