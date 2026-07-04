%% febio
% Below is a demonstration of the features of the |febio| function

%%
clear; close all; clc;

%% Syntax
% |febio(varargin);|

%% Description 
% This function calls FEBio from MATLAB and runs it in the command window.
% If no input is provided this commend simply triggers FEBio to run without
% an input file. 
% If the input is a valid FEBio model input file path then FEBio will run
% the analysis for that file. 
% The input may also be a structure with the following default fields: 
%   defaultOptionStruct.FEBioPath=getFEBioPath; %The FEBio path
%   defaultOptionStruct.run_filename=[]; %The .feb file name
%   defaultOptionStruct.run_logname=[]; %The log file name
%   defaultOptionStruct.runMode='internal'; %Running internally

%% Examples 
% 
% |febio|

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
