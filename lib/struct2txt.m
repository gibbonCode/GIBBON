function struct2txt(varargin)

%% parse input

switch nargin
    case 2
        S=varargin{1};
        fileName=varargin{2};
        delimiter=[];
    case 3
        S=varargin{1};
        fileName=varargin{2};
        delimiter=varargin{3};
end

if isempty(delimiter)
    delimiter = ',';
end

%%

fileID = fopen(fileName,'w');

nameCell=fieldnames(S);
for q=1:1:numel(nameCell)
    entryNow=S.(nameCell{q});
    if isnumeric(entryNow)
        entryText=vec2strIntDouble(entryNow,'%6.7e');
    else
        entryText=entryNow;
    end
    textLine=[nameCell{q},delimiter,entryText];
    fprintf(fileID,[textLine,'\n']);
end
fclose(fileID);

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
