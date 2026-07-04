%% isapprox
% Below is a demonstration of the features of the |isapprox| function

%%
clear; close all; clc;

%% Syntax
% |[L]=isapprox(A,B,tolLevel);|

%% Description 
% This function returns if A is approximately equal to B (to within
% tolLevel) in the form of the logical L. 
% In other words the function simply evaluates: L=abs(A-B)<tolLevel;

%% Examples 
% 

tolLevel=1e-4; % The tolerance level to use

A = - pi; %Example number
B = A + tolLevel/100; % Small difference
C = A + tolLevel*100; % Big difference 
D = A + tolLevel; % Difference equal to the tolerance level

%%
% Using isapprox here should return a true since the difference is smaller
% than the tolerance.
L1 = isapprox(A,B,tolLevel)

%%
% Using isapprox here should return a false since the difference is larger
% than the tolerance.
L2 = isapprox(A,C,tolLevel)

%%
% Using isapprox here should return a false since the difference is equal
% not smaller than the tolerance.
L3 = isapprox(A,D,tolLevel)

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
