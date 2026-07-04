%% cellPair2struct
% Below is a demonstration of the features of the |cellPair2struct| function

%%
clear; close all; clc;

%% Syntax
% |[S]=cellPair2struct(fieldNameCell,fieldDataCell,convertOption,optionStruct);|

%% Description 
% Converts a cell array pair to a structure. The cell |fieldNameCell|
% contains field names and the cell |fieldDataCell| contains the data. The
% data may be anything a structure can hold. The option parameter
% |convertOption| can be used to convert numeric data to text. The optional
% |optionStruct| structure can be used to control this conversion (see
% |mat2strIntDouble|). 

%% Examples 
% 

fieldNameCell={'bob','mary','banana','qwerty_yuiop'};
fieldDataCell={5,'sure',[pi 1 3.5],[1 2 3]};

%%
% Convert all numerical data
convertOption=1;
[S]=cellPair2struct(fieldNameCell,fieldDataCell,convertOption)

%%
% No conversion
convertOption=0;
[S]=cellPair2struct(fieldNameCell,fieldDataCell,convertOption)

%%
% Convert some of the entries
convertOption=[1 0 0 1];
[S]=cellPair2struct(fieldNameCell,fieldDataCell,convertOption)


%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
