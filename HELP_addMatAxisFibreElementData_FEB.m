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