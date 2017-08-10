function [T]=txtfile2cell(fileName)
% function [T]=txtfile2cell(fileName)
% ------------------------------------------------------------------------
% This function read the text file specified by fileName (path and name)
% whereby each line is read into a seperate entry in the output cell array
% T. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/102/10
%------------------------------------------------------------------------

fid=fopen(fileName);
T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
T=T{1,1};
fclose(fid);
 
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
