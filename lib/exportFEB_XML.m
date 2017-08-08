function exportFEB_XML(save_name,XDOC)


% %%
% %Write to text file using xmlwrite (adds extra lines for unknown reasons)
% xmlwrite(save_name,XDOC);
% 
% %Import back into cell array
% [T]=txtfile2cell(save_name);
% 
% %Now save to txt file while skipping "empty lines"
% cell2txtfile(save_name,T,1);

write_XML_no_extra_lines(save_name,XDOC);



 
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
