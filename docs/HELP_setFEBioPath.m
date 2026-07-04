%% setFEBioPath
% Below is a demonstration of the features of the |setFEBioPath| function

%%
clear; close all; clc;

%% Syntax
% |setFEBioPath(FEBioPathSpec);|

%% Description 
% This function sets the FEBio path defined by FEBioPathSpec. Setting the
% path means it is stored in the GIBBON configuration file:
% /GIBBON/config/FEBioPath.txt. 
%
% See also |getFEBioPath|

%% Examples 
% This example shows how the febio path can be specified and set for use in
% GIBBON

%% 
% The FEBio path should be the full path to the FEBio executable file (and
% should include the extendion for Windows). 

%%
% These are examples of paths defined for Windows and Linux (Ubuntu):
% Windows: 
% C:\Program Files\febio-2.9.1\bin\FEBio2.exe
%
% Linux: 
% /home/kevin/FEBio-2.9.1/bin/febio2

%%
% To avoid this demo from undoing the path settings for a user the current
% path is here instead uploaded with |getFEBioPath|

%Define FEBio path
febioPath=getFEBioPath; % Example: febioPath='/home/kevin/FEBio-2.9.1/bin/febio2'

%%
% Setting the FEBio path
setFEBioPath(febioPath); 

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
