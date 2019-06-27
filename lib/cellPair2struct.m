function [S]=cellPair2struct(varargin)

% function [S]=cellPair2struct(fieldNameCell,fieldDataCell,convertOption,optionStruct)
% ------------------------------------------------------------------------
% Converts a cell array pair to a structure. The cell |fieldNameCell|
% contains field names and the cell |fieldDataCell| contains the data. The
% data may be anything a structure can hold. The option parameter
% |convertOption| can be used to convert numeric data to text. The optional
% |optionStruct| structure can be used to control this conversion (see
% |mat2strIntDouble|). 
%
% 2019/06/27 Updated input handling
% 2019/06/27 Added extra conversion options based on mat2strIntDouble
% 2019/06/27 Updated description
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        fieldNameCell=varargin{1};
        fieldDataCell=varargin{2};
        convertOption=zeros(1,numel(fieldNameCell));
        optionStruct=[];
    case 3
        fieldNameCell=varargin{1}; 
        fieldDataCell=varargin{2};
        convertOption=varargin{3}; 
        optionStruct=[];
    case 4
        fieldNameCell=varargin{1};
        fieldDataCell=varargin{2};
        convertOption=varargin{3};
        optionStruct=varargin{4};
end

if numel(convertOption)==1
    convertOption=convertOption*ones(1,numel(fieldNameCell));
end

defaultOptionStruct.formatDouble='%6.7e';
defaultOptionStruct.formatInteger='%d';
defaultOptionStruct.dlmChar=',';
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%%
% Convert cell to structure
S=[];
for q=1:numel(fieldNameCell)
    if convertOption(q)==1 && isnumeric(fieldDataCell{q})
        S.(fieldNameCell{q})=mat2strIntDouble(fieldDataCell{q},optionStruct);
    else
        S.(fieldNameCell{q})=fieldDataCell{q};
    end
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
