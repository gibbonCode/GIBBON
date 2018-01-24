%% febioStruct2xml
% Below is a demonstration of the features of the |febioStruct2xml| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=febioStruct2xml(febio_spec,fileName,optionStruct);|

%% Description
% This function provides the basis for coding FEBio XML input file content
% and helps to generate the XML input file. 

%%
% See also: |febioStructTemplate|

%% About FEBio XML input files
% FEBio input files (.feb) are XML files. The XML syntax bellow is for
% febio_spec 2.5. More on the input file formatting can be found in the
% FEBio user manual.

%% Coding XML files in MATLAB
% Coding XML files (or in this case .feb FEBio files) can be challenging
% and may look confusing. An alternative method is implemented here which
% is a trade-off between clarify, ease of programming, and fidelity with
% the generated XML file. In general XML files follow the following form:
 
%%

% <Globals> 
%   <Constants>
%       <R>8.314e-6</R>
%       <T>298</T>
%       <Fc>96485e-9</Fc>
%   </Constants>
%   <Solutes>
%       <solute id="1" name="Na">
%           <charge_number>1</charge_number>
%           <molar_mass>22.99</molar_mass>
%           <density>0.97<density>
%       </solute>
% </Globals>

%%
% Here Globals is an element with the entries Constants and Solutes being
% child elements (and Constants and Solutes have R,T,Fc and solute as their
% respective child elements,... etc). The symbols <..> start an element and
% </..> closes that element. Components within <..> sometimes contain
% attributes. An example of this is id="1", here the attribute "id" has the
% value "1" associated with it, another example of an attribute is
% name="Na". Elements may contain child elements and may have a value. For
% instance the element Constants has children but the element R has the
% value 8.314e-6 associated with it. 
% To code these features in MATLAB a particular syntax was developed.
% Structures are suitable as they are tree-like similar to XML. However
% structures lack the attribute feature, i.e. fields in structures have
% values but do not have attributes associated with them. To solve this
% keywords are used here to denote attributes and to distinguish them from
% values or elements. For instance node.ATTR.id=1 will become <node
% id="1">, where the ATTR keyword is here simply used as a target/trigger
% keyword for determining wheter the subsequent feature is an element,
% value or attribute. Futhermore the VAL keyword exists to denote a value,
% e.g. node.ATTR.id=1; node.VAL=2;  will become:

% <node id="1"> 2 </node>. 

%%
% One may use arrays for values as well i.e. node.ATTR.id=nodeIdArray;
% node.VAL=nodeCoordinateArray; can be used to define sets of attributes
% and sets of values simultaneously. 
% Below several types of FEBio functionality and how to code for the XML
% content is highlighted. Examples are obtained or adjusted from the FEBio
% user manual. Not all of FEBio's functionality is covered in this
% documentation example. However in principle all FEBio functionality
% (based on XML content) can be coded for in a similar fashion. 

%% Set febio_spec version
febio_spec.ATTR.version='2.5';

%% Module section
%%
% XML syntax example
%%

% <Module type="Solid"/>
%%
% Implementation example
febio_spec.Module.ATTR.type='solid'; %Use default set

%% Control section

%Example structure
defaultStruct.analysis.type='static';
defaultStruct.title='temp';
defaultStruct.time_steps=10;
defaultStruct.step_size=0.1;
defaultStruct.dtol=0.001;
defaultStruct.etol=0.01;
defaultStruct.rtol=0;
defaultStruct.lstol=0.9;
defaultStruct.time_stepper.dtmin=defaultStruct.step_size/3;
defaultStruct.time_stepper.dtmax=defaultStruct.step_size*3; %Constant
% defaultStruct.time_stepper.dtmax.ATTR.lc=1; %Load curve based
% defaultStruct.time_stepper.dtmax.VAL=0.1; %Load curve based
defaultStruct.time_stepper.max_retries=5;
defaultStruct.time_stepper.opt_iter=10;
defaultStruct.max_refs=15;
defaultStruct.max_ups=10;
defaultStruct.optimize_bw=0;
% defaultStruct.restart.ATTR.file='restartDumpFile.dmp'; %On with file specified
defaultStruct.restart.VAL=0; %Off
defaultStruct.plot_level='PLOT_MAJOR_ITRS';
defaultStruct.plot_range=[0,-1];
defaultStruct.plot_stride=1;
defaultStruct.plot_zero_state=0;
defaultStruct.cmax=1e5;
defaultStruct.print_level='PRINT_MINOR_ITRS';
defaultStruct.min_residual=1e-20;
% defaultStruct.integration='N/A';
defaultStruct.output_level='OUTPUT_MAJOR_ITRS'; 

febio_spec.Control=defaultStruct; %Use default set

%% Globals section

% <Globals> 
%   <Constants>
%       <R>8.314e-6</R>
%       <T>298</T>
%       <Fc>96485e-9</Fc>
%   </Constants>
%   <Solutes>
%       <solute id="1" name="Na">
%           <charge_number>1</charge_number>
%           <molar_mass>22.99</molar_mass>
%           <density>0.97<density>
%       </solute>
%       <solute id="2" name="Cl">
%           <charge_number>-1</charge_number>
%           <molar_mass>35.45</molar_mass>
%           <density>3.21<density>
%       </solute>
%       <solute id="3" name="Glc">
%           <charge_number>0</charge_number>
%           <molar_mass>180.16</molar_mass>
%           <density>1.54<density>
%       </solute>
%   </Solutes>
% </Globals>

febio_spec.Globals.Constants.T=0; 
febio_spec.Globals.Constants.R=0;
febio_spec.Globals.Constants.Fc=0; 

febio_spec.Globals.Solutes{1}.ATTR.id=1; 
febio_spec.Globals.Solutes{1}.ATTR.name='Na';
febio_spec.Globals.Solutes{1}.charge_number=1;
febio_spec.Globals.Solutes{1}.molar_mass=22.99;
febio_spec.Globals.Solutes{1}.density=0.97;

febio_spec.Globals.Solutes{2}.ATTR.id=2; 
febio_spec.Globals.Solutes{2}.ATTR.name='Cl';
febio_spec.Globals.Solutes{2}.charge_number=-1;
febio_spec.Globals.Solutes{2}.molar_mass=35.45;
febio_spec.Globals.Solutes{2}.density=3.21;

%% Material section

febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=1;
febio_spec.Material.material{1}.m1=6;
febio_spec.Material.material{1}.c2=1;
febio_spec.Material.material{1}.m2=-6;
febio_spec.Material.material{1}.k=100;

%% Geometry section

%%
% Nodes

% <Nodes name="set01">
%   <node id="1">0,0,0,</node>
%   ...
%   </node id="101">1,1,1</node>
%   </Nodes>
%   <Nodes name ="set02">
%   <node id="102">2,1,1</node>
%   ...
%   <node id="999">2,2,2</node>
% </Nodes>

%Nodes which are defined with a set name
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet1';
febio_spec.Geometry.Nodes{1}.node.ATTR.id=[1 102 2]';
febio_spec.Geometry.Nodes{1}.node.VAL=randn(3,3);

%Nodes which are defined with a set name
febio_spec.Geometry.Nodes{2}.ATTR.name='nodeSet2';
febio_spec.Geometry.Nodes{2}.node.ATTR.id=(12:18)';
febio_spec.Geometry.Nodes{2}.node.VAL=rand(7,3);

%Nodes without a set name
% febio_spec.Geometry.Nodes{3}.ATTR.name='nodeSet3';
n=7; %Number of nodes to test
febio_spec.Geometry.Nodes{3}.node.ATTR.id=(1:n)';
febio_spec.Geometry.Nodes{3}.node.VAL=rand(n,3);

%%
% Node sets

% <NodeSet name="nodeset1">
%     <node id="1"/>
%     <node id="102"/>
%     <node id="2"/>
% </NodeSet>
febio_spec.Geometry.NodeSet{1}.ATTR.name='nodeSet3';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=(1:21)';

% <NodeSet name="set1"> ... </NodeSet>
febio_spec.Geometry.NodeSet{2}.ATTR.name='nodeSet4';
febio_spec.Geometry.NodeSet{2}.VAL=(1:21);

%%
% Elements

% <Elements type="hex8" mat="1" name="part1">
%   <elem id="1">1,2,3,4,5,6,7,8</elem>
%   <elem id="2">1,2,3,4,5,6,7,8</elem>
%   <...>
% </Elements>

febio_spec.Geometry.Elements{1}.ATTR.type='hex8';
febio_spec.Geometry.Elements{1}.ATTR.mat=1;
febio_spec.Geometry.Elements{1}.ATTR.name='part1';
febio_spec.Geometry.Elements{1}.elem.ATTR.id=2+(1:12)';
febio_spec.Geometry.Elements{1}.elem.VAL=round(100*randn(12,8));

%%
% Edges

% <Edge name="edge1">
%   <line2 lid="1">1,2</line2>
%   <line2 lid="2">2,3</line2>
% </Edge>

febio_spec.Geometry.Edge{1}.ATTR.name='edge1';
febio_spec.Geometry.Edge{1}.line2.ATTR.lid=(1:2)';
febio_spec.Geometry.Edge{1}.line2.VAL=round(100*rand(2,2));

febio_spec.Geometry.Edge{2}.ATTR.name='edge1';
febio_spec.Geometry.Edge{2}.line3.ATTR.lid=(1:3)';
febio_spec.Geometry.Edge{2}.line3.VAL=round(100*rand(3,2));

%%
% Surfaces

% %Surface section
% % <Surface name="named_surface">
% % <quad4 lid="1">1,2,3,4</quad4>
% % <...>
% % </Surface>

febio_spec.Geometry.Surface{1}.ATTR.name='Surface1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:4)';
febio_spec.Geometry.Surface{1}.quad4.VAL=round(100*rand(4,4));

febio_spec.Geometry.Surface{2}.ATTR.name='Surface2';
febio_spec.Geometry.Surface{2}.tri7.ATTR.lid=(1:5)';
febio_spec.Geometry.Surface{2}.tri7.VAL=round(100*rand(5,7));

%%
% Discrete

% <DiscreteSet name="springs">
%   <delem>1,2</delem>
%   <delem>3,4</delem>
% </DiscreteSet>

febio_spec.Geometry.DiscreteSet{1}.ATTR.name='Springs';
febio_spec.Geometry.DiscreteSet{1}.delem.VAL=round(100*rand(3,2));

%%
% ElementSet

% <ElementSet name="set01">
%   <elem id="1001"/>
%   <elem id="1002"/>
%   <elem id="1003"/>
% </ElementSet>

febio_spec.Geometry.ElementSet{1}.ATTR.name='ElementSet1';
febio_spec.Geometry.ElementSet{1}.elem.ATTR.id=(1:5)';

%%
% SurfacePair

% <SurfacePair name="contact1">
%   <master surface="surface1"/>
%   <slave surface="surface2"/>
% </SurfacePair>

febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface='Suface1';
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface='Surface2';

%% MeshData section

%%

% <MeshData>
%   <ElementData var="shell thickness" elem_set="part1">
%       <elem lid="1">1.0, 1.0, 1.0, 1.0</elem>
%       <elem lid="2">1.0, 1.0, 1.0, 1.0</elem>
%       </ElementData >
%   <ElementData var="pre_stretch" elem_set="part2">
%       <elem lid="1">1.05</elem>
%       <elem lid="2">1.05</elem>
%   </ElementData >
% </MeshData>

febio_spec.MeshData.ElementData{1}.ATTR.var='Shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set='part1';
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:2)';
febio_spec.MeshData.ElementData{1}.elem.VAL=round(100*rand(2,7));

febio_spec.MeshData.ElementData{2}.ATTR.var='pre_stretch';
febio_spec.MeshData.ElementData{2}.ATTR.elem_set='part2';
febio_spec.MeshData.ElementData{2}.elem.ATTR.lid=(1:2)';
febio_spec.MeshData.ElementData{2}.elem.VAL=[1.05 1.05]';

%% Initial section

%%

% <MeshData>
%   <NodeData name="init_vel" node_set="set1">
%       <node lid="1">1.0</node>
%       <node lid="2">1.5</node>
%       <node lid="3">1.0</node>
%   </NodeData>
% </MeshData>
% <Initial>
%   <init bc="vx" node_set="set1">
%       <value node_data="init_vel"/>
%   </init>
% </Initial>

febio_spec.MeshData.NodeData{1}.ATTR.name='init_vel';
febio_spec.MeshData.NodeData{1}.ATTR.node_set='set1';
febio_spec.MeshData.NodeData{1}.node.ATTR.lid=(1:3)';
febio_spec.MeshData.NodeData{1}.node.VAL=[1.0 1.5 1.0]';

febio_spec.Initial.init.ATTR.bc='vx';
febio_spec.Initial.init.ATTR.node_set='set1';
febio_spec.Initial.init.value.ATTR.node_data='init_vel';

%% Boundary

% <MeshData>
%   <NodeData name="values" node_set="set1">
%       <node lid="1">1.0</node>
%       <node lid="2">2.0</node>
%       <node lid="3">3.0</node>
%   </NodeData>
% </MeshData>
%
% <Boundary>
%   <prescribe bc="x" node_set="set1">
%       <scale lc="1">2.0</scale>
%       <relative>0</relative>
%       <value node_data="values"/>
%   </prescribe>
%   <fix bc="x" node_set="nodeset1"/>
%   <rigid rb="2" node_set="set1"/>
% </Boundary>

febio_spec.MeshData.NodeData{2}.ATTR.name='init_vel';
febio_spec.MeshData.NodeData{2}.ATTR.node_set='set1';
febio_spec.MeshData.NodeData{2}.node.ATTR.lid=(1:3)';
febio_spec.MeshData.NodeData{2}.node.VAL=[1.0 1.5 1.0]';

febio_spec.Boundary.prescribe{1}.ATTR.bc='x';
febio_spec.Boundary.prescribe{1}.ATTR.node_set='set1';
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=2;
febio_spec.Boundary.prescribe{1}.relative=0;
febio_spec.Boundary.prescribe{1}.value.ATTR.node_data='values';

febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set='nodeset1';

febio_spec.Boundary.rigid{1}.ATTR.rb=2;
febio_spec.Boundary.rigid{1}.ATTR.node_set='nodeset1';

%%
% Prescribed Rigid Body Degrees of Freedom

% <rigid_body mat="1">
%   <prescribed bc="x" lc"1">2.0</prescribed>
%   <fixed bc="y"/>
%   <fixed bc="Rx"/>
% </rigid_body>

febio_spec.Boundary.rigid_body{1}.ATTR.mat=1;
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed.VAL=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='Rx';

%% Loads section 

%%
% Nodal loads

% <Loads>
%   <nodal_load bc="x" node_set="set1">
%       <scale lc="1">1.0</scale>
%       <value>3.14</value>
%   </nodal_load>
% </Loads>

febio_spec.Loads.nodal_load.ATTR.bc='x';
febio_spec.Loads.nodal_load.ATTR.node_set='set1';
febio_spec.Loads.nodal_load.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load.scale.VAL=1;
febio_spec.Loads.nodal_load.value.VAL=pi;

%% LoadData section

% <LoadData>
%   <loadcurve id="1" [type="type" extend="extend"]>
%       <point> 0, 0 </point>
%       <point> 1, 1 </point>
%   </loadcurve>
% </LoadData>

febio_spec.LoadData.loadcurve.ATTR.id=1;
febio_spec.LoadData.loadcurve.ATTR.type='linear';
febio_spec.LoadData.loadcurve.point.VAL=[0 0; 1 1];

%% Output section 

%%
% Logfile

% <Output>
%   <logfile [file="<log file>"]>
%       <node_data [attributes]>item list</node_data>
%       <element_data [attributes]>item list</element_data>
%       <rigid_body_data [attributes]>item list</rigid_body_data>
%   </logfile>
% </Output>
febio_spec.Output.logfile.ATTR.file='logfile.txt';

febio_spec.Output.logfile.node_data{1}.ATTR.file='outputLogfile1.txt';
febio_spec.Output.logfile.node_data{1}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:10;

febio_spec.Output.logfile.element_data{1}.ATTR.file='outputLogfile1.txt';
febio_spec.Output.logfile.element_data{1}.ATTR.data='sxx;syy;szz;sxy;syz;sxz';
febio_spec.Output.logfile.element_data{1}.ATTR.name='Element stresses';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:25;

febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file='outputLogfile3.txt';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx;Fy;Fz';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';

%%
% Plot file

% <Output>
%   <plotfile type="febio">
%       <var type="displacement"/>
%       <var type="stress"/>
%   </plotfile>
% </Output>

febio_spec.Output.plotfile.ATTR.type='febio';
febio_spec.Output.plotfile.var{1}.ATTR.type='displacement';
febio_spec.Output.plotfile.var{2}.ATTR.type='stress';

%% Creating the FEBio input file
% You can use |febioStruct2xml| to write the xml data to a file AND/OR a
% domNode object (optional output). Leave the fileName variable empty to
% supress file export. 

%Create file name for XML file
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
fileName=fullfile(savePath,'tempModel.feb');

[domNode]=febioStruct2xml(febio_spec,fileName); %Exporting to file and domNode

%% Viewing the FEBio input file
% The |febView| command can be used to render and XML file in a figure
% window. 
%%
% NOTE: The figure below does not render in documentation due to a
% MATLAB but (or limitation). The code |[hFig]=febView(domNode);| is
% therefore suppressed. (see also |xmlView|);

% [hFig]=febView(fileName);

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

