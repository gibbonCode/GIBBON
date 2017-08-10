function cleanUpTetGen(pathNameTempFiles)

% function cleanUpTetGen(pathNameTempFiles)
% ------------------------------------------------------------------------
% This function can clean a folder from tetgen output files. 
% 
% See also |cleanDir|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------

extCell={'ele','node','face','edge','mtr','smesh','p2t'}; %Extensions of files to delete

%Remove the files with matching extensions
cleanDir(pathNameTempFiles,extCell);

 
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
