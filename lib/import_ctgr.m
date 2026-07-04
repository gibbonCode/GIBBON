function [outputStruct]=import_ctgr(fileName)

% [outputStruct]=import_ctgr(ctgrFileName)
% ------------------------------------------------------------------------
%
%
% Change log: 
% 2021/10/10 Created
%
%------------------------------------------------------------------------

%% Access XML object
ctgrXML = xmlread(fileName);

%% Parse XML object

timestep_level = ctgrXML.getElementsByTagName('timestep');

numTimeSteps=timestep_level.getLength; %Number of time steps

for iTime=1:1:numTimeSteps    
    contour_level = timestep_level.item(iTime-1).getElementsByTagName('contour');
    numContours=contour_level.getLength;
    for iContour=1:1:numContours
        contour_points_level = contour_level.item(iContour-1).getElementsByTagName('contour_points');
        point_level = contour_points_level.item(0).getElementsByTagName('point');
        num_point=point_level.getLength;
        
        V_points=zeros(num_point,3); %Allocate array for current point set        
        for iPoint=1:1:num_point
            % id_points(iPoint)=str2double(point_level.item(iPoint-1).getAttribute('id').toCharArray()');
            V_points(iPoint,:)=[sscanf(point_level.item(iPoint-1).getAttribute('x').toCharArray()','%f')...
                                sscanf(point_level.item(iPoint-1).getAttribute('y').toCharArray()','%f')...
                                sscanf(point_level.item(iPoint-1).getAttribute('z').toCharArray()','%f')];
        end        
        outputStruct.contour{iContour}.contour_points=V_points;        
    end
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
