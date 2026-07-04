%% DCM2euler
% Below is a demonstration of the features of the |DCM2euler| function

%%
clear; close all; clc;

%% Syntax
% |[a]=DCM2euler(Q);|

%% Description 
% This function is the inverse of the euler2DCM function. The Euler angles
% |a| are derived based on the input rotatin tensor |Q|. 

%% Examples 
% 

%%
% Plot settings
fontSize=25;

%% Retrieving the Euler andles from a rotation tensor

%% 
% Get example patch data
[F,V]=parasaurolophus;

%%
% Defining sets of true Euler angles for X, Y and Z axis rotation
a_true=[0.25*pi 0.25*pi -0.25*pi]

%%
% Use |euler2DCM| function to define the rotation tensor
[Q]=euler2DCM(a_true);

%%
% Use |DCM2euler| to retrieve the Euler angles
a_fit=DCM2euler(Q)

%% Handling symbolic expressions

try
    syms a b c
    
    a_true=[a b c]
    Q=euler2DCM(a_true);
    
    a_fit=DCM2euler(Q)
    
catch
    warning('Symbolic toolbox likely missing')
end

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
