function export_off(fileName,F,V)


%%
%Get edges
[E]=patchEdges(F,1);

%Create top text
T={'OFF ';...
       [num2str(size(V,1)),' ',num2str(size(F,1)),' ',num2str(size(E,1))];...
       ''};
   
%Create vertex text
textForm=repmat('%f ',1,size(V,2));
textForm=[textForm,'\n'];
TV=sprintf(textForm,V');
T{end+1}=TV(1:end-1); %Take off added end of line statement

%Create faces text
F_mat=[size(F,2)*ones(size(F,1),1) F-1]; 
textForm=repmat('%u ',1,size(F,2)+1);
textForm=[textForm,'\n'];
TF=sprintf(textForm,F_mat');
T{end+1}=TF;

%Write text to file
cell2txtfile(fileName,T,0);   
 
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
