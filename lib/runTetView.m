function [varargout]=runTetView(modelName)

%% SETTING TETGEN PATHNAMES

compString=computer; 
switch compString
    case 'PCWIN' %Windows 32-bit
        error('PCWIN 32-bit is not supported. Compile tetGen from the source and alter the code here');
    case 'PCWIN64' %Windows 64-bit
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','win64');
        runNameTetView=fullfile(pathNameTetView,'tetview-win.exe');
    case 'GLNXA64'        
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','lin64');
        runNameTetView=fullfile(pathNameTetView,'tetview-linux');
    case 'MACI64'        
        error('MACI64 is not supported yet. Get TetView online and alter the code here');
    otherwise
        error('Your platform does not seem to be supported. Code your own solution or contact support.')
end

modelName=regexprep(modelName,'\','/');

%% RUN TETVIEW

runString=['"',runNameTetView,'" "',modelName,'" & '];
[runStatus,runCmdHist]=system(runString);

varargout{1}=runStatus;
varargout{2}=runCmdHist;

 
%% 
% ********** _license boilerplate_ **********
% 
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
