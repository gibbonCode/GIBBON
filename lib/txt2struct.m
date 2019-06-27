function [S]=txt2struct(varargin)

%% parse input

switch nargin
    case 1
        fileName=varargin{1};
        delimiter=[];
    case 2
        fileName=varargin{1};
        delimiter=varargin{2};
end

if isempty(delimiter)
    delimiter = ',';
end
%%
if exist(fileName,'file')
    %Import data into cell
    
    formatSpec = '%s%[^\n\r]';
    fileID = fopen(fileName,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    % Convert cell to structure        
    nameCell=dataArray{1};
    entryCell=dataArray{2};    
    [S]=cellPair2struct(nameCell,entryCell,1);
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
