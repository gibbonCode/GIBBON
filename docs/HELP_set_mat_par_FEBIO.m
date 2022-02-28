%% set_mat_par_FEBIO
% Below is a demonstration of the features of the |set_mat_par_FEBIO| function

%%
clear; close all; clc;

%%
% Folder/file locations
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','FEB');

%% Altering material parameter entries
% First we specify the feb file either by its file name or the XML object,
% the latter is shown here

febFileName=fullfile(pathName,'tempModel_2p0.feb');

%% 
% View .FEB file using the |xmlView| function
xmlView(febFileName);

%%
% Lets assume the material section looks like this: 

%    <Material>
%       <material id="1" name="mat_1" type="Ogden">
%          <c1>1.0000000e-03</c1>
%          <m1>1.2000000e+01</m1>
%          <k>1.0000000e+00</k>
%       </material>
%       <material id="2" name="mat_2" type="rigid body">
%          <density>1.0000000e+00</density>
%          <center_of_mass>0.0000000e+00, 0.0000000e+00, 0.0000000e+00</center_of_mass>
%       </material>
%    </Material>

%%
% But this is desired (c1, k and density are altered):

%    <Material>
%       <material id="1" name="mat_1" type="Ogden">
%          <c1>3.1415927e+00</c1>
%          <m1>1.2000000e+01</m1>
%          <k>3.1415927e+03</k>
%       </material>
%       <material id="2" name="mat_2" type="rigid body">
%          <density>1.7320508e+00</density>
%          <center_of_mass>0.0000000e+00, 0.0000000e+00, 0.0000000e+00</center_of_mass>
%       </material>
%    </Material>

%%
% Then a cell array is composed containig an entry per material that
% requires altered parameters. Each cell centry contains a structure with
% the fields id, par_names and par_values. For the example above one would
% define:

%Entry for material id 1
mat_struct.id=1;
mat_struct.par_names={'c1','k'};
mat_struct.par_values={pi,pi*1e3};
mat_cell{1}=mat_struct;

%Entry for material id 2
mat_struct.id=2;
mat_struct.par_names={'density'};
mat_struct.par_values={1.2345};
mat_cell{2}=mat_struct;

%% 
% Then the function |set_mat_par_FEBIO| is called to replace the material
% parameter entries

febFileSaveName=fullfile(pathName,'tempModel.feb'); %The new file name

set_mat_par_FEBIO(febFileName,febFileSaveName,mat_cell);

%%
% Verify the XML file content
xmlView(febFileSaveName);

%% Altering material parameter entries for "nested" or multisolid materials

febFileName=fullfile(pathName,'tempModel_2p0_solidMixture.feb');

xmlView(febFileName);

%%
% Lets assume the material section looks like this: 

%    <Material>
%      <material id="1" name="mat_1" type="solid mixture">
%          <solid type="Ogden">
%             <c1>8.4582737e-04</c1>
%             <m1>1.1907214e+01</m1>
%             <k>1.1734536e+00</k>
%          </solid>
%          <solid type="fiber-exp-pow-uncoupled">
%             <mat_axis type="user"/>
%             <ksi>1.5010797e-03</ksi>
%             <alpha>1.0000000e-25</alpha>
%             <beta>3.0527617e+00</beta>
%             <theta>0.0000000e+00</theta>
%             <phi>0.0000000e+00</phi>
%          </solid>
%       </material>
%       <material id="2" name="mat_2" type="rigid body">
%          <density>1.7320508e+00</density>
%          <center_of_mass>0.0000000e+00, 0.0000000e+00, 0.0000000e+00</center_of_mass>
%       </material>
%    </Material>

%%
% But this is desired (c1, ksi and density are altered):

%    <Material>
%       <material id="1" name="mat_1" type="solid mixture">
%          <solid type="Ogden">
%             <c1>3.1415927e+00</c1>
%             <m1>1.1907214e+01</m1>
%             <k>1.1734536e+00</k>
%          </solid>
%          <solid type="fiber-exp-pow-uncoupled">
%             <mat_axis type="user"/>
%             <ksi>1.2345000e+00</ksi>
%             <alpha>1.0000000e-25</alpha>
%             <beta>3.0527617e+00</beta>
%             <theta>0.0000000e+00</theta>
%             <phi>0.0000000e+00</phi>
%          </solid>
%       </material>
%       <material id="2" name="mat_2" type="rigid body">
%          <density>4.5678000e+00</density>
%          <center_of_mass>0.0000000e+00, 0.0000000e+00, 0.0000000e+00</center_of_mass>
%       </material>
%    </Material>

%%
% This time another cell array is composed, again containig an entry per
% material that requires altered parameters. However this time the
% par_names field is different since now it contains cell entries
% containing first references to the solid followed by the actual parameter
% name. For the example above one would define:

%Entry for material id 1
mat_struct.id=1;
mat_struct.par_names={{'solid','Ogden','c1'},{'solid','fiber-exp-pow-uncoupled','ksi'}};
mat_struct.par_values={pi,1.2345}; 
mat_cell{1}=mat_struct;

%Entry for material id 2
mat_struct.id=2;
mat_struct.par_names={'density'};
mat_struct.par_values={4.5678};
mat_cell{2}=mat_struct;

%% 
% Then the function |set_mat_par_FEBIO| is called to replace the material
% parameter entries

febFileSaveName=fullfile(pathName,'tempModel.feb'); %The new file name

set_mat_par_FEBIO(febFileName,febFileSaveName,mat_cell);

%%
% Verify the XML file content
xmlView(febFileSaveName);

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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
