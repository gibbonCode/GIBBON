function export_INP(elementStruct,nodeStruct,saveName)

% export_INP(elementStruct,nodeStruct,saveName)
% ------------------------------------------------------------------------
% CURRENTLY THIS FUNCTION ONLY SUPPORTS 1 ELEMENT TYPE PER EXPORT NO MIXED
% ELEMENT MODELS!
%
% This function exports ABAQUS INP input files using the following inputs:
% 
% nodeStruct : A structure array containing the following fields: 
%       nodeStruct.N :An nx3 array of nodal coordinates
%       nodeStruct.N_ind :An nx1 array of the node indices (numbers)
%
% elementStruct : A a structure array of the form (for example):
%       E_type: '*ELEMENT, TYPE=STRI3, ELSET=PART-DEFAULT_1_EB1'  
%               Specifying element type as the ABAQUS string
%       E: [4536x3 double] An mxl array of the nodal connectivity
%       E_ind: [4536x1 double] An mx1 array for the element indices
%
% EXAMPLE switch structure for element type  
% switch elementStruct.E_type
%     case 'tri3'
%         elementStruct.E_type='*ELEMENT, TYPE=STRI3, ELSET=PART-DEFAULT_1_EB1';
%     case 'quad4'
%         elementStruct.E_type='*ELEMENT, TYPE=S4R, ELSET=PART-DEFAULT_1_EB1';
%     case 'tet4'
%         elementStruct.E_type='*ELEMENT, TYPE=C3D4, ELSET=PART-DEFAULT_1_EB1';
%     case 'hex8'
%         elementStruct.E_type='*ELEMENT, TYPE=C3D8R, ELSET=PART-DEFAULT_1_EB1';
% end
%
% saveName : the file name for the new INP file. The file extension should
% be .inp. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 23/08/2012
%------------------------------------------------------------------------

%%
disp('--- export_INP ---');

%% LOAD DEFAULT FILE STRUCTURE

%Loading mother INP file
[T]=txtfile2cell('motherINP.inp'); %This can be any INP file containing the targets

% Finding node and element fields in mother INP file
targets={'*NODE','*ELEMENT','**'}; %Text form targets
lineCount=1;
lineIndexTarget=zeros(size(targets));
i_target=1;
while 1
    l=T{lineCount};
    target=targets{i_target};
    if (strfind(l,target))
        lineIndexTarget(i_target)=lineCount;
        i_target=i_target+1;
    end    
    lineCount=lineCount+1;
    if nnz(lineIndexTarget)==(numel(targets))
        break
    end
end

%% CREATING NODE AND ELEMENT TEXT FIELDS
disp('CREATING NODE AND ELEMENT TEXT FIELDS');

disp('---> Creating node text field');
%Creating new node text field
NODE_field=cell(size(nodeStruct.N,1),1);
for q=1:1:size(nodeStruct.N,1)
    nodeFieldLine=[sprintf('%8d,',nodeStruct.N_ind(q)) sprintf('% 10.6e, ',nodeStruct.N(q,:)) ];   
    NODE_field{q}=nodeFieldLine(1:end-2);    
end
 
%Creating new element text field
disp('---> Creating element text field');
ELEMENT_field=cell(size(elementStruct.E,1),1);
for q=1:1:size(elementStruct.E,1)    
    elementFieldLine=sprintf('%8d,',[elementStruct.E_ind(q) elementStruct.E(q,:)]);   
    ELEMENT_field{q}=elementFieldLine(1:end-1);    
end

%% CREATING NEW TEXT CELL

%Change element type line
T(lineIndexTarget(2))={elementStruct.E_type};

%Update node and element text fields
T=[T(1:lineIndexTarget(1)); NODE_field; T(lineIndexTarget(2)-2:lineIndexTarget(2)); ELEMENT_field; T(lineIndexTarget(3):end)];

%% WRITTING INP FILE
disp('EXPORTING TO INP FILE...');

fid=fopen(saveName,'wt');
for q=1:size(T,1)
    fprintf(fid,'%s\n',T{q});
end
fclose(fid);

disp('DONE!');
end
 
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
