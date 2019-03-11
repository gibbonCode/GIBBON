function [varargout]=importNodeFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %f %f %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
nodeID=double(A{1});
V=nan(max(nodeID),3);
V(nodeID,:)=[A{2} A{3} A{4}];

varargout{1}=nodeID;
varargout{2}=V;

 
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
