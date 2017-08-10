function [pathNames]=getSubPaths(pathName)

% function [pathNames]=getSubPaths(pathName)
% ------------------------------------------------------------------------
%
%This function (based on the GENPATH command) creates the output pathNames
%which contains all folders and sub-folders within the folder specified by
%the input pathName.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/04/2013
% 2017/06/01 Fixed bug in relation to operational system differences
%------------------------------------------------------------------------

%%
if ispc % semi-colon splits paths
    strPattern=[filesep,';'];    
else %colon splits paths for Linux and mac
    strPattern=':';        
end
pathNames=regexp(genpath(pathName),strPattern, 'split');
pathNames=pathNames(2:end-1)';

end
 
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
