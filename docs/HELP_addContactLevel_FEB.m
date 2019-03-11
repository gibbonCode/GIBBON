%% addContactLevel_FEB
% Below is a demonstration of the features of the |addContactLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addContactLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the contact information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Defining contact
% 

%Example data 
contactPenalty=500; 
contactType=1; 
F1=[1 2 3]; %Faces for surface 1
F2=[4 5 6]; %Faces for surface 2

%Defining surfaces
FEB_struct.Geometry.Surface{1}.Set=F1;
FEB_struct.Geometry.Surface{1}.Type='tri3';
FEB_struct.Geometry.Surface{1}.Name='Contact_master_indentor';

FEB_struct.Geometry.Surface{2}.Set=F2;
FEB_struct.Geometry.Surface{2}.Type='tri3';
FEB_struct.Geometry.Surface{2}.Name='Contact_slave_gel';

%Adding contact information
FEB_struct.Contact{1}.Surface{1}.SetName=FEB_struct.Geometry.Surface{1}.Name;
FEB_struct.Contact{1}.Surface{1}.Type='master';
FEB_struct.Contact{1}.Surface{2}.SetName=FEB_struct.Geometry.Surface{2}.Name;
FEB_struct.Contact{1}.Surface{2}.Type='slave';
switch contactType
    case 1 %STICKY
        FEB_struct.Contact{1}.Type='sticky';
        FEB_struct.Contact{1}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance','max_traction','snap_tol'};
        FEB_struct.Contact{1}.Values={0,0.1,contactPenalty,0,10,0.01,0,0.01};
    case 2 %SLIDING facet-facet
        FEB_struct.Contact{1}.Type='facet-to-facet sliding';
        FEB_struct.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Contact{1}.Values={contactPenalty,1,0,...
            0,0.1,...
            0,0,10,...
            0,...
            0.01};
    case 3 %SLIDING sliding_with_gaps
        FEB_struct.Contact{1}.Type='sliding_with_gaps';
        FEB_struct.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'fric_coeff','fric_penalty',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Contact{1}.Values={contactPenalty,1,0,...
            0,0.1,...
            0,0,10,...
            0,1,...
            0,...
            0.01};
end

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

% %Add geometry information (Surfaces)
% domNode=addGeometryLevel_FEB(domNode,FEB_struct);

%Add boundary condition information
domNode=addContactLevel_FEB(domNode,FEB_struct);

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
