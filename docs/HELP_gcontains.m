%% gcontains
% Below is a demonstration of the features of the |gcontains| function

%%
clear; close all; clc;

%% Syntax
% |TF=contains(str,strPattern);|

%% Description
% This function is an alternative to the MATLAB function |contains|,
% introduced in MATLAB R2016b. The function attempts to use MATLAB
% |contains| but uses custom implementation if |contains| is not found. See
% also: |contains|.

%% Examples

%% String array

str = ["Mary Ann Jones","Christopher Matthew Burns","John Paul Smith"];

pattern = ["Ann","Paul"];
% TF = contains(str,pattern)
TF = gcontains(str,pattern)

%% Using |IgnoreCase|

str = ["Anne","Elizabeth","Marianne","Tracy"];

pattern = "anne";
% TF = contains(str,pattern,'IgnoreCase',true)

TF = gcontains(str,pattern,'IgnoreCase',true)

%% Character arrays

chr = 'peppers, onions, and mushrooms';

% TF = contains(chr,'onion')
TF = gcontains(chr,'onion')

% TF = contains(chr,'pineapples')
TF = gcontains(chr,'pineapples')
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
