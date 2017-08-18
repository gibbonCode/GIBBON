%% addConstraintsLevel_FEB
% Below is a demonstration of the features of the |addConstraintsLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addConstraintsLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the contraints information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Defining fix/prescribed constraint for a rigid body
% 

%Example data 
bcPrescribeMagnitudes=[pi 1 -5]; %example prescribed displacement magnitudes
rigidMaterialId=2; %Example material Id for rigid body material

%Constraint section
FEB_struct.Constraints{1}.RigidId=rigidMaterialId;

FEB_struct.Constraints{1}.Prescribe{1}.bc='x';
FEB_struct.Constraints{1}.Prescribe{1}.Scale=bcPrescribeMagnitudes(1);
FEB_struct.Constraints{1}.Prescribe{1}.lc=1;

FEB_struct.Constraints{1}.Prescribe{2}.bc='y';
FEB_struct.Constraints{1}.Prescribe{2}.Scale=bcPrescribeMagnitudes(2);
FEB_struct.Constraints{1}.Prescribe{2}.lc=1;

FEB_struct.Constraints{1}.Prescribe{3}.bc='z';
FEB_struct.Constraints{1}.Prescribe{3}.Scale=bcPrescribeMagnitudes(3);
FEB_struct.Constraints{1}.Prescribe{3}.lc=1;

FEB_struct.Constraints{1}.Fix{1}.bc='Rx';
FEB_struct.Constraints{1}.Fix{2}.bc='Ry';
FEB_struct.Constraints{1}.Fix{3}.bc='Rz';

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addConstraintsLevel_FEB(domNode,FEB_struct);

%%
%  View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
