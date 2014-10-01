function [F,V] = import_off(fileName)
%
% OFF
% numVertices numFaces 0
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
numFacesVertices = str2num(fgets(fid));
numVertices=numFacesVertices(1);
numFaces=numFacesVertices(2);

%Get vertices
[V] = fscanf(fid,'%f %f %f',[3 numVertices])';

% read Face 1  1088 480 1022
[F_set] = fscanf(fid,'%d %d %d %d\n',[4 numFaces])';
F = F_set(:,2:end)+1;

fclose(fid);

