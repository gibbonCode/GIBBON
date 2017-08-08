function [varargout]=importNodeFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %f %f %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
nodeID=double(A{1});
V=nan(max(nodeID),3);
V(nodeID,:)=[A{2} A{3} A{4}];

varargout{1}=nodeID;
varargout{2}=V;

 
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
