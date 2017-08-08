function FEBioPath=getFEBioPath

% function FEBioPath=getFEBioPath
% ------------------------------------------------------------------------
% This function reads a user specified FEBio path name from a config text
% file. The text file may contain multiple paths (e.g. for multiple
% operational systems), however the first valid patch is always used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/10/05 Updated for GIBBON
% 2015/10/05 Enabled searching valid paths in multi-line config. file.
% 2017/06/10 changed empty output to an empty character array. 
%------------------------------------------------------------------------

%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths
try 
    [T]=txtfile2cell(fileName);    
catch
    T={};
end

%Get first valid path
FEBioPath='';
if ~isempty(T)
    for q=1:1:numel(T)
        pathName=T{q};
        if exist(pathName,'file')==2
            FEBioPath=pathName;
            return
        end
    end
end

 
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
