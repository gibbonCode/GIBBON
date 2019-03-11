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
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
