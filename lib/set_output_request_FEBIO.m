function set_output_request_FEBIO(filename,savename,logfile_names,output_types,data_types)

% Adds output section for writting data to logfile. The function does NOT
% check if the output types already exists and thus the new output request
% is always apended. 
%
% Example:
%
% <febio_spec>
% ...................
%    <Output>
%       ...................
%       <logfile>
%           <[output_type]="[data_type]" delim=", " file="[logfile_name]"/>
%       </logfile>
%    </Output>
% ...................
% </febio_spec>
%
% 25/01/2012, Kevin Mattheus Moerman
% kevinmoerman@gmail.com

%%
%Get XML
docNode = xmlread(filename);

%Get docRootNode
FEBIOSPEC_Node = docNode.getElementsByTagName('febio_spec').item(0);
docRootNode =  FEBIOSPEC_Node;

%Check for Output level
if docNode.getElementsByTagName('Output').getLength==0; %Output level does not exist yet
    new_output=1;    
    %Adding Output level
    OutputNode = docNode.createElement('Output');
    docRootNode = docRootNode.appendChild(OutputNode);        
else %Output level already exists
    new_output=0;    
    docRootNode=docNode.getElementsByTagName('Output').item(0);   
end

%Check for existing logfile entries
no_entries=docRootNode.getElementsByTagName('logfile').getLength;
if no_entries==0 %No existing logfile entries
    new_logfile=1;    
    %Adding logfile level
    logfile_node = docNode.createElement('logfile');
    logfile_node = docRootNode.appendChild(logfile_node);    
else %Logfile entries exist
    new_logfile=0;
    logfile_node=docRootNode.getElementsByTagName('logfile').item(0);
end

%Adding output requests
if ~iscell(output_types)
    output_types={output_types};
    logfile_names={logfile_names};
    data_types={data_types};
end
for q=1:1:numel(output_types)
    
    output_type=output_types{q};
    logfile_name=logfile_names{q};
    data_type=data_types{q};
    
    %Adding output_type level
    output_type_node = docNode.createElement(output_type);
    logfile_node.appendChild(output_type_node);
    
    %Adding output_type attribute data_type
    attr = docNode.createAttribute('data');
    attr.setNodeValue(data_type);
    output_type_node.setAttributeNode(attr);
    
    %Adding output_type attribute logfile_name
    attr = docNode.createAttribute('file');
    attr.setNodeValue(logfile_name);
    output_type_node.setAttributeNode(attr);
    
    %Adding output_type attribute delim
    attr = docNode.createAttribute('delim');
    attr.setNodeValue(', ');
    output_type_node.setAttributeNode(attr);
    
end

%% Saving XML file
% xmlwrite(savename,docNode);
write_XML_no_extra_lines(savename,docNode);

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
