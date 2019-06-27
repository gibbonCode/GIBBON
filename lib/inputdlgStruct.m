function S=inputdlgStruct(varargin)

%% Parse input
switch nargin
    case 1        
        nameCell=varargin{1};
        defaultOptions=[];
        dialogTitle=[];
    case 2        
        nameCell=varargin{1};
        defaultOptions=varargin{2};
        dialogTitle=[];
    case 3
        nameCell=varargin{1};
        defaultOptions=varargin{2};
        dialogTitle=varargin{3};
end

if isempty(dialogTitle)
   dialogTitle='Input dialog'; 
end

if isempty(defaultOptions)
   defaultOptions=repmat({''},size(nameCell));
end
  
%%
s=25+max([cellfun(@numel,nameCell) cellfun(@numel,defaultOptions)]); %Set sizes
entryCell = inputdlg(nameCell,dialogTitle,[1 s],defaultOptions); %Open dialog box
[S]=cellPair2struct(nameCell,entryCell,1); %Convert output to structure

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
