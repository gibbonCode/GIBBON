function [varargout]=febStruct2febFile(FEB_struct)

% function [varargout]=febStruct2febFile(FEB_struct)
% ------------------------------------------------------------------------
%
% This function uses the input FEB_struct to generate an FEBio .feb file.
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/05/27: Updated for GIBBON
% 2015/05/09: Added biphasic capabilities
% 2015/05/10: Added traction load capabilities
% 2017/08/18: Adding support to allow incomplete feb file creation (only
% write fields provided). 
%------------------------------------------------------------------------

%%

dispStartTitleGibbonCode('Writing FEBio XML object');

%% Set default display setting if missing
if ~isfield(FEB_struct,'disp_opt')
    FEB_struct.disp_opt=0;
end

%% Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field

%Set febio_spec
febio_spec = domNode.getDocumentElement;
if ~isfield(FEB_struct,'febio_spec')
    FEB_struct.febio_spec.version='2.0';
elseif ~isfield(FEB_struct.febio_spec,'version')
    FEB_struct.febio_spec.version='2.0';
end
febio_spec.setAttribute('version',FEB_struct.febio_spec.version); %Adding version attribute

%% Add comment if present
if isfield(FEB_struct,'commentField')
    commentString = FEB_struct.commentField;
else %Default comment
    commentString = ['Created using GIBBON, ',datestr(now)];
end
commentNode = domNode.createComment(commentString);
febio_spec.appendChild(commentNode);

%% DEFINING MODULE LEVEL
if ~isfield(FEB_struct,'Module')
    FEB_struct.Module.Type='solid'; %Use solid as default module
end
domNode=addModuleLevel_FEB(domNode,FEB_struct);

%% DEFINE CONTROL SECTION
if isfield(FEB_struct,'Control')
    domNode=addControlLevel_FEB(domNode,FEB_struct);
end

%% DEFINING GLOBALS LEVEL
domNode=addGlobalsLevel_FEB(domNode,FEB_struct);

%% DEFINING MATERIAL LEVEL
domNode=addMaterialLevel_FEB(domNode,FEB_struct);

%% DEFINING GEOMETRY LEVEL
writeMethod=1;
if isfield(FEB_struct,'Geometry')
    switch writeMethod
        case 1 % TEXT FILE PARSING (faster for large arrays)
            domNode=addGeometryLevel_TXT(domNode,FEB_struct);
        case 2 %XML PARSING
            domNode=addGeometryLevel_FEB(domNode,FEB_struct);
    end
end

%% DEFINE BOUNDARY CONDITIONS LEVEL AND LOADS LEVEL
if isfield(FEB_struct,'Boundary')
    [domNode]=addBoundaryLevel_FEB(domNode,FEB_struct);
end

%% DEFINE LOAD LEVEL
if isfield(FEB_struct,'Loads')
    [domNode]=addLoadsLevel_FEB(domNode,FEB_struct);
end

%% DEFINE CONTACT
if isfield(FEB_struct,'Contact')
    [domNode]=addContactLevel_FEB(domNode,FEB_struct);
end

%% DEFINE CONSTRAINTS LEVEL
if isfield(FEB_struct,'Constraints')
    [domNode]=addConstraintsLevel_FEB(domNode,FEB_struct);
end

%% DEFINE LOADDATA LEVEL
if isfield(FEB_struct,'LoadData')
    [domNode]=addLoadDataLevel_FEB(domNode,FEB_struct);
end

%% DEFINE STEP LEVEL
if isfield(FEB_struct,'Step')
    [domNode]=addStepLevel_FEB(domNode,FEB_struct);
end

%% DEFINE OUTPUT LEVEL
domNode=addOutputLevel_FEB(domNode,FEB_struct);

%% CREATE OUTPUT OR EXPORT XML FILE

switch nargout
    case 0
        disp('Writing .feb file');
        if isfield(FEB_struct,'topCommentLine')
            exportFEB_XML(FEB_struct.run_filename,domNode,FEB_struct.topCommentLine); % Saving XML file
        else
            exportFEB_XML(FEB_struct.run_filename,domNode); % Saving XML file
        end
    case 1
        varargout{1}=domNode;
end
dispDoneGibbonCode;

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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
