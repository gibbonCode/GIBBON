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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
