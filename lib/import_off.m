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
