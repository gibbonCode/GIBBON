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
%
testGibbon('help','test',1,2);

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
