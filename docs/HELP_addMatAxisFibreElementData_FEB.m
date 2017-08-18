%% addMatAxisFibreElementData_FEB
% Below is a demonstration of the features of the |addMatAxisFibreElementData_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addMatAxisFibreElementData_FEB(domNode,FEB_struct);|

%% Description
% This function adds the local element matetrial axis or fiber direction
% information to the input XML data (domNode) based on the input febio
% structure (FEB_struct).

%% Examples

%% Example: Defining the local material axis for each element
% 

n=5; %number of elements
E=round(20*rand(n,8)); %Simulates element matrix
V_fib=zeros(n,3); %Simulated fibre directions
V_fib(:,3)=1; 

%Adding fibre direction, construct local orthonormal basis vectors
[a,d]=vectorOrthogonalPair(V_fib);

VF_E=zeros(size(V_fib,1),size(V_fib,2),2);
VF_E(:,:,1)=a; %a1 ~ e1 ~ X or first direction
VF_E(:,:,2)=d; %a2 ~ e2 ~ Y or second direction
%Vf_E %a3 ~ e3 ~ Z, third direction, e.g. fibre direction

FEB_struct.Geometry.Nodes=rand(25,3);
FEB_struct.Geometry.Elements={E}; %The element sets
FEB_struct.Geometry.ElementType={'hex8'}; %The element types
FEB_struct.Geometry.ElementMat={ones(size(E,1),1)};
FEB_struct.Geometry.ElementsPartName={'Block'};

FEB_struct.Geometry.ElementData.MatAxis.ElementIndices=1:1:size(V_fib,1);
FEB_struct.Geometry.ElementData.MatAxis.Basis=VF_E;

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add geometry information (NodeSet)
domNode=addGeometryLevel_FEB(domNode,FEB_struct);
   
%Add boundary condition information
domNode=addMatAxisFibreElementData_FEB(domNode,FEB_struct);

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
