function [F,V] = import_off(fileName)

% function [F,V] = import_off(fileName)
% ------------------------------------------------------------------------
% This function imports http://www.geomview.org/ type .off files specifying
% a surface mesh (provided all faces are of the same type!). The output
% consists of (patch data) faces F and vertices V. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------

%%

%Get file ID
fid = fopen(fileName,'r');

%Check if file is valid
if( fid==-1 )
    error('File not found!');
end

%Check first line of OFF to denote a .off file
firstLine = fgets(fid);  %First line
if isempty(strfind(firstLine, 'OFF'))
    error('First line should contain the OFF keyword to denote a .off file!');    
end

%Access number of faces and vertices
done=0; 
while done==0
    currentLine=fgets(fid);
    if ~isempty(currentLine)
        numFacesVertices = str2num(currentLine);
        if numel(numFacesVertices)==3
            numVertices=numFacesVertices(1);
            numFaces=numFacesVertices(2);
            done=1; 
        end
    end
end

%Get vertices
[V] = fscanf(fid,'%f %f %f\n',[3 numVertices])';

%Read next line which is for the face set
currentPosition = ftell(fid); %Get current position
currentLine=fgets(fid); %First face line
fseek(fid,currentPosition,-1); %Put back position
currentLineDouble=str2num(currentLine); %Current line converted to numbers
numVerticesPerFace=currentLineDouble(1); %First number is number of vertices per face e.g. 3 for a triangle

%Get faces (currently assuming all faces are of the same type!)
t_form=repmat('%d ',1,numVerticesPerFace+1); 
t_form=[t_form(1:end-1),'\n'];
[F_set] = fscanf(fid,t_form,[numVerticesPerFace+1 numFaces])';
F = F_set(:,2:end)+1;

fclose(fid);

 
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
