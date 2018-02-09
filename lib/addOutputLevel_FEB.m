function docNode=addOutputLevel_FEB(docNode,FEB_struct)

disp('Adding Output level')

rootNode = docNode.getDocumentElement;

outputNode = docNode.createElement('Output');
outputNode = rootNode.appendChild(outputNode);

%% Adding plotfile field
disp('----> Adding plotfile field')

parent_node = docNode.createElement('plotfile');
parent_node = outputNode.appendChild(parent_node);
parent_node.setAttribute('type','febio');

if ~isfield(FEB_struct,'Output')
    FEB_struct.Output.VarTypes={'displacement','stress','relative volume'}; %DEFAULT
end

% Adding Output requests
if isfield(FEB_struct,'Output');
    for q=1:1:numel(FEB_struct.Output.VarTypes)
        var_node = docNode.createElement('var'); %create entry
        var_node = parent_node.appendChild(var_node); %add entry
        var_node.setAttribute('type',FEB_struct.Output.VarTypes{q}); %add attribute
    end
end
%% Adding logfile field

if isfield(FEB_struct,'run_output_names')
    
    disp('----> Adding logfile field')
    
    %Adding logfile level
    logfile_node = docNode.createElement('logfile');
    logfile_node = outputNode.appendChild(logfile_node);
    
    % FEB_struct.run_filename,FEB_struct.run_filename,FEB_struct.run_output_names,FEB_struct.output_types,FEB_struct.data_types
    % (filename,savename,logfile_names,output_types,data_types)
    
    %Adding output requests
    output_types=FEB_struct.output_types;
    logfile_names=FEB_struct.run_output_names;
    data_types=FEB_struct.data_types;
    
    if ~iscell(output_types)
        output_types={output_types};
        logfile_names={logfile_names};
        data_types={data_types};
    end
    
    for q=1:1:numel(output_types)
        
        output_type=output_types{q};
        logfile_name=logfile_names{q};
        data_type=data_types{q};

        %Remove path from filename is present
        [pathName,fileName,fileExtension] = fileparts(logfile_name);
        if ~isempty(pathName)
            
            %TEMPORARY FIX TO COPE WITH OLDER FEBIO VERSIONS
            if ~isfield(FEB_struct,'ignoreLogPathIssue')
                removePath=1;
            else
                removePath=~FEB_struct.ignoreLogPathIssue;
            end
            if removePath==1
                warning('Provided path of logfile is replaced by .feb file path. Only provide filename to avoid this warning');
                logfile_name=[fileName,fileExtension]; %Combine only name and extension
            end
            
        end
        
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
