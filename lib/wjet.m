function cMap=wjet(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[98  161 177;...
    130 189 164;...
    157 209 147;...
    186 219 151;...
    211 228 152;...
    231 236 154;...
    240 236 146;...
    249 228 141;...
    255 210 130;...
    251 180 113;...
    248 149  95;...
    243 124  87;...
    231  97  83;...
    216  67  80;...
    ];

cMap=cMap./255;

[cMap]=resampleColormap(cMap,n);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
