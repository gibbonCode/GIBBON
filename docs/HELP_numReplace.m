%% numReplace
% Below is a demonstration of the features of the |numReplace| function

%%
clear; close all; clc;

%% REPLACING NUMBERS IN ARRAYS
%%
% An example array
A=[0,-6,3,0,0;-1,-4,-7,-9,-9;-4,-7,11,-5,12;10,-7,5,-7,13;0,11,-2,-6,2;12,4,2,-5,NaN]

%%
% Defining the input array for entries (NaN's allows) that need to be replaced
a=[2 -5 nan 0]
%%
% Defining the numbers (NaN's allows) to take their place
b=[991 992 993 994]; %Numbers to take their place

%%
% Replacing the numbers using |numReplace|
B=numReplace(A,a,b)

%% NOTES ON PERFORMANCE FOR NON-INTEGERS
% The |numReplace| function employs the |ismember| function. Hence it is
% suitable for all number cases where ismember is able to detect
% membership. Numerical precission difficulties may arise for non-integer
% entires. Consider the below: 

logicMember=ismember(pi,pi+eps(pi))
logicMember=ismember(pi,pi+eps(pi)/10)

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
