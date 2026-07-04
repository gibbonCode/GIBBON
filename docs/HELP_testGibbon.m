%% testGibbon
% Below is a demonstration of the features of the |testGibbon| function

%%
clear; close all; clc;

%% Syntax
% |testGibbon(varargin);|

%% Description 
% |testGibbon| can be used to run the help and demo (documentation) file
% for GIBBON. By default the function will run all help and demo files in
% test mode. However the user may select just the help files or just the
% demo files and choose publish mode instead of test mode. 
%
% Optional inputs:
%        testSet --- 'all', 'help', or 'demo'
%        testMode ---'test' or 'pub' i.e. test run the file or publish the
%        file
%        approveQuestion --- 1 or 0, For 1 the user will be asked to proceed to the next file. If no is answered the file is opened for editing 
%        startInd --- e.g. 1 for the first file 
% 

%% Examples 
testSet='all'; % Or 'demo' or 'help'
testMode='test'; % for test mode or 'pub' for publishing mode
approveQuestion=1; % Ask user to approve the test 1 or more on 0
startInd=1; % File index to start testing from

%%
% |testGibbon(testSet,testMode,approveQuestion,startInd);|

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
