function D=import_FBG(loadName)
%-------------------------------------------------------------------------
%
%function D=import_FBG(loadName)
%
%A simple function to import FBG data from text file (specified by
%loadName) and output the data fields (columns) in a structure array D.
%
%
%Kevin Mattheus Moerman
%kevinmoerman@hotmail.com
%2014/01/07 - Created function
%-------------------------------------------------------------------------

%% PARSE TEXT FILE

fid=fopen(loadName);
[cell_out] = textscan(fid,'%s %f %f %f %f %f %f\n', 'delimiter', ',');
fclose(fid);

%% FORMULATE OUTPUT

D.test_date=cell_out{1};
D.cycle_no=cell_out{2};
D.FBG_temp=cell_out{3};
D.FBG_strain=cell_out{4};
D.ind_speed=cell_out{5};
D.ind_depth=cell_out{6};
D.TTL_logic=cell_out{7};

 
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
