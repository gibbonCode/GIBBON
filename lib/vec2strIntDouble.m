function [t]=vec2strIntDouble(varargin)

%%
switch nargin
    case 1
        n=varargin{1};
        formatDouble='%6.7e';
        rowWrapLength=[];
    case 2
        n=varargin{1};
        formatDouble=varargin{2};
        rowWrapLength=[];
    case 3
        n=varargin{1};
        formatDouble=varargin{2};
        rowWrapLength=varargin{3};
end

%%

optionStruct.formatDouble=formatDouble;
optionStruct.rowWrapLength=rowWrapLength;
[t]=mat2strIntDouble(n,optionStruct);

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
