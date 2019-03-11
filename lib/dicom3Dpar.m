function [varargout]=dicom3Dpar(D)

%Gets the 3D parameters from the DICOM info file
try
    v=D.PixelSpacing;
catch
    v=[1 1];
    warning('Missing field PixelSpacing, assuming [1 1]!');
end

%N.B. TO DO TREAT SPACING/THICKNESS PROPERLY (currently suitable for 3D
%volumes) 
try
    v(3)=D.SpacingBetweenSlices;
catch
    try
        v(3)=D.SliceThickness;
    catch
        v(3)=1;
        warning('Missing field SpacingBetweenSlices and SliceThickness, assuming 1!');
    end
end

% ORIGIN:
% The x, y, and z coordinates of the upper left hand corner (center of the
% first voxel transmitted) of the image, in mm. 
try
    OR=D.ImagePositionPatient;
catch
    OR=[0 0 0];
    warning('Missing field ImagePositionPatient, assuming [0 0 0]!');
end

try
    r=D.ImageOrientationPatient(4:6); % ROW DIRECTION
    c=D.ImageOrientationPatient(1:3); % COLUMN DIRECTION
catch
    r=[0 1 0];
    c=[1 0 0];
    warning('Missing field ImageOrientationPatient, assuming [1 0 0 0 1 0]!');
end

switch nargout
    case 1 %If only one output argument is requested then create structure
        G.v=v;
        G.OR=OR;
        G.r=r;
        G.c=c;
        varargout{1}=G; 
    otherwise
        varargout{1}=v;
        varargout{2}=OR;
        varargout{3}=r;
        varargout{4}=c;
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
