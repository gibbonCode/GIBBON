%% prin2effective
% Below is a demonstration of the features of the |prin2effective| function

%%
clear; close all; clc;

%% Syntax
% |[D_eff]=prin2effective(D1,D2,D3,typeFlag);|
% |[D_eff]=prin2effective(D,typeFlag);|

%% Description 
% This function computes the effective/Von Mises stress or strain based on
% the input principal components. If typeFlag='stress' then the von Mises
% stress is computer. If typeflag='strain' the effective strain is
% computed. 

%% Examples 
% 

%% Example 1: Computing von Mises stress from a triplet of principal stresses

% A zero Von Mises stress state
S_prin=[2 2 2]; 
S_vm=prin2effective(S_prin,'stress')

% A unit Von Mises stress state
S_prin=[1 2 1]; 
S_vm=prin2effective(S_prin,'stress')

% A sqrt(3) Von Mises stress state
S_prin=[1 2 3]; 
S_vm=prin2effective(S_prin,'stress')

%% Example 2: Computing effective strain from a triplet of principal strains

% A zero effective strain state
E_prin=[0.5 0.5 0.5]; 
E_effective=prin2effective(E_prin,'strain')

% A unit effective strain state
E_prin=[0.5 2 2]; 
E_effective=prin2effective(E_prin,'strain')

% A 1/3 effective strain state
E_prin=[-0.25 0.25 0.25]; 
E_effective=prin2effective(E_prin,'strain')

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
