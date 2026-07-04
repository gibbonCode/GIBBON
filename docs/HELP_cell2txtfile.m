%% cell2txtfile
% Below is a demonstration of the features of the |cell2txtfile| function

%% Syntax
% |cell2txtfile(fileName,T,skipOpt);|

%% Description
% The |cell2txtfile| function exports the content in the cell array T to
% the text file fileName. Each entry in the cell array will be a line in
% the txt file. Prior to text file creation the cell is converted to a
% column format. If the input skipOpt=1 cell entries which appear empty
% (after spaces are removed) will be skipped.

%%
clear; close all; clc;

%% Examples

%% Exporting a cell containing text data to a txt file  
% Create example cell containing text entries, empty entries, an entry
% containing just spaces and a non-char entry number. The number entry will
% be converted to a string. 

T={'Hello','','World','  ',uint8(125),pi};
filePath=mfilename('fullpath');
fileName=fullfile(fileparts(fileparts(filePath)),'data','temp','temp.txt');
skipOpt=0; %Empty entries will be kept 
cell2txtfile(fileName,T,skipOpt);

%%
%Output text reread as cell:
[T_out]=txtfile2cell(fileName);
disp(T_out)

%%
% Changing skipOpt to 1 will skip these empty lines
skipOpt=1; %Empty entries will be skipped
cell2txtfile(fileName,T,skipOpt);

[T_out]=txtfile2cell(fileName);
disp(T_out)

%%
% Remove example file
delete(fileName);

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
