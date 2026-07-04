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
% e.g. node.ATTR.id=1; node.VAL=2;  will become: <node id="1"> 2 </node>. 

%%
% One may use arrays for values as well i.e. node.ATTR.id=nodeIdArray;
% node.VAL=nodeCoordinateArray; can be used to define sets of attributes
% and sets of values simultaneously. 
%% 
% In principle all FEBio functionality (based on XML content) can be coded
% for in a similar fashion. Study the FEBio user manual for more
% information. 

%%
% The remainder of this documentation will describe how various XML
% components of the FEBio input file can be coded. A structure is created
% with the name "febio_spec". For more information on coding the febio
% structure see the various FEBio demos provided. 

%% Set febio_spec version
febio_spec.ATTR.version='4.0';

%% Module section
% XML syntax example:
% <Module type="Solid"/>

%%
% Implementation example
febio_spec.Module.ATTR.type='solid'; 

%% Control section


%% Example 1: Globals section

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

%% Example 2: Material section

%%
% An uncouple material
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=1;
febio_spec.Material.material{1}.m1=6;
febio_spec.Material.material{1}.k=100;

%%
% A visco-elastic material

%Viscoelastic part
febio_spec.Material.material{2}.ATTR.type='uncoupled viscoelastic';
febio_spec.Material.material{2}.ATTR.Name='Block_material';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.g1=0.1;
febio_spec.Material.material{2}.t1=0.01;
febio_spec.Material.material{2}.density=1e-9;

%Elastic part
febio_spec.Material.material{2}.elastic{1}.ATTR.type='Ogden';
febio_spec.Material.material{2}.elastic{1}.c1=1;
febio_spec.Material.material{2}.elastic{1}.m1=2;
febio_spec.Material.material{2}.elastic{1}.k=100;
febio_spec.Material.material{2}.elastic{1}.density=1e-9;

%%
% A poro-elastic material

%Viscous part
febio_spec.Material.material{3}.ATTR.type='biphasic';
febio_spec.Material.material{3}.ATTR.name='Block_material';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.phi0=0.5;
febio_spec.Material.material{3}.permeability.ATTR.type='perm-const-iso';
febio_spec.Material.material{3}.permeability.ATTR.name='permeability';
febio_spec.Material.material{3}.permeability.perm=7.41e-11; %hydraulic permeability
febio_spec.Material.material{3}.fluid_density=1e-9;

%Solid part
febio_spec.Material.material{3}.solid{1}.ATTR.type='Ogden unconstrained';
febio_spec.Material.material{3}.solid{1}.c1=1;
febio_spec.Material.material{3}.solid{1}.m1=2;
febio_spec.Material.material{3}.solid{1}.cp=100;
febio_spec.Material.material{3}.solid{1}.density=1e-9;

%%
% A poro-elastic and visco-elastic material

% 	<Material>
% 		<material id="1" name="Biphasic-Viscoelastic" type="biphasic">
% 			<phi0>0.2</phi0>
% 			<solid type="viscoelastic">
% 				<g1>0.5</g1>
% 				<t1>3000</t1>
% 				<elastic type="solid mixture">
% 					<solid type="neo-Hookean">
% 						<density>1</density>
% 						<E>0.5</E>
% 						<v>0</v>
% 					</solid>
% 					<solid type="ellipsoidal fiber distribution">
% 						<beta>2,2,2</beta>
% 						<ksi>5,5,5</ksi>
% 					</solid>
% 				</elastic>
% 			</solid>
% 			<permeability type="perm-Holmes-Mow">
% 				<perm>0.001</perm>
% 				<M>3</M>
% 				<alpha>2</alpha>
% 			</permeability>
% 		</material>

%Viscous part
febio_spec.Material.material{4}.ATTR.type='biphasic';
febio_spec.Material.material{4}.ATTR.name='Biphasic-Viscoelastic';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.phi0=0.2;
febio_spec.Material.material{4}.permeability.ATTR.type='perm-Holmes-Mow';
febio_spec.Material.material{4}.permeability.perm=0.001; 
febio_spec.Material.material{4}.permeability.M=3; 
febio_spec.Material.material{4}.permeability.alpha=2; 

%Solid part
%-> Viscoelastic part
febio_spec.Material.material{4}.solid{1}.ATTR.type='viscoelastic';
febio_spec.Material.material{4}.solid{1}.g1=0.5;
febio_spec.Material.material{4}.solid{1}.t1=3000;
febio_spec.Material.material{4}.solid{1}.density=1e-9;

%-> Elastic part
febio_spec.Material.material{4}.solid{1}.elastic{1}.ATTR.type='solid mixture';
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{1}.density=1e-9;
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{1}.E=0.5;
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{1}.v=0;
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{2}.ATTR.type='ellipsoidal fiber distribution';
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{2}.beta=[2 2 2];
febio_spec.Material.material{4}.solid{1}.elastic{1}.solid{2}.ksi=[5 5 5];
         
%% Example 3: Mesh section

%Example geometry
[meshStruct]=hexMeshBox([2 2 2],[2 2 2]);
E=meshStruct.E; V=meshStruct.V;
bcSupportList=(1:size(V,1))';

% -> Nodes

% <Nodes name="Object1">
%   <node id="1">0,0,0,</node>
%   <node id="2">1,0,0,</node>
%   ...
% </Nodes>

febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements

% <Elements type="hex8" name="Part1">
%   <elem id="1">1,2,3,4,5,6,7,8</elem>
%   <elem id="2">1,2,3,4,9,10,11,12</elem>
%   <...>
% </Elements>

partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%% Creating the FEBio input file
% You can use |febioStruct2xml| to write the xml data to a file AND/OR a
% domNode object (optional output). Leave the fileName variable empty to
% supress file export. 

%Create file name for XML file
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
fileName=fullfile(savePath,'tempModel.feb');

optionStruct.arrayParseMethod=1;
[domNode]=febioStruct2xml(febio_spec,fileName,optionStruct);

%% Viewing the FEBio input file
% The |febView| command can be used to render an XML file in a figure
% window. Alternatively the |textView| command can be used:
h=febView(fileName,1);

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
