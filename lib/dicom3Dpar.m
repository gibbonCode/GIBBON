function [v,OR,r,c]=dicom3Dpar(D)

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

if nargout==1 %If only one output argument is requested then create structure
   G.v=v; 
   G.OR=OR;
   G.r=r;
   G.c=c;
   v=G;
end

