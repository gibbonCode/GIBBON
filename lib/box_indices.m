function [IND]=box_indices(siz)

% function [IND]=box_indices(siz)
% ------------------------------------------------------------------------
% This function returns the indices of the boundary entries (e.g. pixels or
% voxels) for an array with size siz. The boundary entries form a box at
% the edge of the matrix. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

numDim=numel(siz); 
logicSides=false(siz); 
permuteOrder=1:numDim;
for q=1:1:numDim
    if q>1
        logicSides=permute(logicSides,permuteOrder);
    end
    
    switch numDim
        case 2
            logicSides(1,:)=1;
            logicSides(end,:)=1;
        case 3
            logicSides(1,:,:)=1;
            logicSides(end,:,:)=1;
        case 4
            logicSides(1,:,:,:)=1;
            logicSides(end,:,:,:)=1;
    end
       
    if q>1
        logicSides=ipermute(logicSides,permuteOrder);
    end
    
    permuteOrder=[permuteOrder(2:end) permuteOrder(1)]; %Shift order
   
end

IND=find(logicSides);
 
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
