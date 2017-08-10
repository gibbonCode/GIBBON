function write_XML_no_extra_lines(save_name,XDOC)

% XML_string = xmlwrite(XDOC); %XML as string
% 
% XML_string = regexprep(XML_string,'\n[ \t\n]*\n','\n'); %removes extra tabs, spaces and extra lines
% 
% %Write to file
% fid = fopen(save_name,'w');
% fprintf(fid,'%s\n',XML_string);
% fclose(fid); 

%Write to text file
%xmlwrite(save_name,XDOC); %MATLAB's xmlwrite extremely slow for v2017a-2017b
xmlwrite_xerces(save_name,XDOC); %Custom XML write function for now

%Import back into cell array
[T]=txtfile2cell(save_name);

%Save to txt file while skipping "empty lines"
cell2txtfile(save_name,T,1);
 
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
