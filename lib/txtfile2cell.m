function [T]=txtfile2cell(fileName)
% function [T]=txtfile2cell(fileName)
% ------------------------------------------------------------------------
% This function read the text file specified by fileName (path and name)
% whereby each line is read into a seperate entry in the output cell array
% T. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/102/10
%------------------------------------------------------------------------

fid=fopen(fileName);
T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
T=T{1,1};
fclose(fid);
 
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
