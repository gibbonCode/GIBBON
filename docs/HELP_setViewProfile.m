%% setViewProfile
% Below is a demonstration of the features of the |setViewProfile| function

%%
clear; close all; clc;

%% Syntax
% |setViewProfile(profileName);|

%% Description 
% This function sets the vcw view manipulation profile to use. The
% following are currently supported: 
% * CAD
% * febio
% * touchpad
%
% See also: |getViewProfile|, |vcw|

%% Examples 
% 

%% Example 1: Setting a view profile
% The options 'CAD', 'febio', and 'touchpad' are available. 

% Using currently set value as example here to avoid this demo from
% overwriting a user set value
profileName=getViewProfile %Desired profile

%%
% Setting the profile
setViewProfile(profileName);

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
