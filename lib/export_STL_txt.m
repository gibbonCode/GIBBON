function export_STL_txt(fileName,stlStruct)

% ------------------------------------------------------------------------
% function export_STL_txt(fileName,stlStruct)
%
% This function generates an STL file using the input structure stlStruct.
% Multiple solids are supported as cell entries within the structure.
%
% Example if the structure array for a single solid:
%
% stlStruct =
%    solidNames: {'stanford_bunny'} %Cell containing solid names
% solidVertices: {[8745x3 double]} %Cell constaining vertices for each solid
%    solidFaces: {[2915x3 double]} %Cell containing faces for each solid
%  solidNormals: {[2915x3 double]} %Cell containing normals for each solid
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% CHANGE LOG:
% 2014/04/08 Created
% 2016/01/13 Added check for file deletion
% ------------------------------------------------------------------------

%% Write STL file in txt format

%Remove file if it exists (in some cases data is appended to file if it
%exists already)
if exist(fileName,'file')==2
    delete(fileName);    
    %Check if its gone
    if exist(fileName,'file')==2
        error(['Existing file with name: ',fileName,' found. Deletion not succesful, check user/file permissions']);
    end
end

%Append STLs to file
writeOpt='a';
for q=1:1:numel(stlStruct.solidNames)   
    patch2STL(fileName,stlStruct.solidVertices{q},stlStruct.solidFaces{q},stlStruct.solidNormals{q},stlStruct.solidNames{q},writeOpt);
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
