function [varargout]=importFaceFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
faceID=double(A{1});
F=nan(max(faceID),3);
F(faceID,:)=double([A{2} A{3} A{4}]);
faceBoundaryID=double(A{5});

if all(isnan(faceBoundaryID(:)))
    faceBoundaryID=-ones(size(F,1),1);
end

varargout{1}=faceID;
varargout{2}=F;
varargout{3}=faceBoundaryID;

 
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
