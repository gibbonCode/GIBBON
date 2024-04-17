%% quiverTriad
% Below is a demonstration of the features of the |quiverTriad| function

%%
clear; close all; clc;

%% Syntax
% |[h]=quiverTriad(V,R,vecSize,colorOpt,alphaLevel);|

%% Description

%% Examples

R = eye(3,3);
V = [0 0 0];

%%
vecSize = 10;

%%
figure;
h(1)=quiverTriad(V,R,vecSize,'r',1);
h(2)=quiverTriad(V,R,vecSize,'r',1);
axisGeom;
drawnow;



%
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
