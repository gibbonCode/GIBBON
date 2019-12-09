function patch2STL(varargin)

% ------------------------------------------------------------------------
% function patch2STL(fileName,V,F,N,solidName)
% 
% This function generates an STL file using the vertices (V) and triangular
% faces (F) provided
% 
% 
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/08 
%
% CHANGE LOG: 
% 09/12/2009
% * Created
% 2014/04/08 
% * Changed %f to %6.6e
% * Added solidName input option
% * Created varargin based input parsing
% ------------------------------------------------------------------------

%% Parse input

%Check if at least 3 inputs and not more than 6 are provided
if nargin<3 || nargin>6
    error('Wrong number of input arguments');
end

%Get first 3 manditory inputs
fileName=varargin{1};
V=varargin{2};
F=varargin{3};

% Check if faces input is consistent with triangles
if ~size(F,2)==3
    error('Patch data does not represent triangles. Convert to triangles first.');
end

%Check N and solidName entries
switch nargin
    case 3 %N, solidName and writeOpt not provided
        N=[];
        solidName='path2STL';
        writeOpt='w';
    case 4 %solidName and writeOpt not provided
        N=varargin{4};        
        solidName='path2STL';
        writeOpt='w';
    case 5
        N=varargin{4};
        solidName=varargin{5};        
        writeOpt='w';
    case 6
        N=varargin{4};
        solidName=varargin{5};        
        writeOpt=varargin{6};        
end

%Check if the solidName is a string
if ~ischar(solidName)
    error('solidName needs to be a string');
end

%% Get face normals
if isempty(N) %if face normals are missing
    [N,~]=trinorm(F(:,[1 2 3]),V);
end

%% Write STL file in ASCII format

fid=fopen(fileName,writeOpt); %Append to existing file

fprintf(fid,['solid ',solidName,' \n']);
hw = waitbar(0,'Writing STL-file, please wait...');

%Start loop for all faces, vertices and normals
for qFace = 1:1:size(F,1);
    
    %Write face normals
    fprintf(fid,'  facet normal %6.6e %6.6e %6.6e\n',N(qFace,:));
    fprintf(fid,'    outer loop\n');
    fprintf(fid,'      vertex %6.6e %6.6e %6.6e\n',V(F(qFace,1),:));
    fprintf(fid,'      vertex %6.6e %6.6e %6.6e\n',V(F(qFace,2),:));
    fprintf(fid,'      vertex %6.6e %6.6e %6.6e\n',V(F(qFace,3),:));
    fprintf(fid,'    endloop\n');
    fprintf(fid,'  endfacet\n');
    
    %Increment waitbar
    waitbar(qFace/size(F,1),hw,['Writing STL-file ',num2str(round((qFace/size(F,1)).*100)),'%']); 
end
fprintf(fid,['endsolid ',solidName,' \n']);
waitbar(1,hw,'Finished! Closing file.');
drawnow;
fclose(fid); %close file ID
close(hw); %close waitbar

 
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
