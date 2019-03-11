function [varargout]=febStruct2febFile(FEBioInputStruct)

% function [varargout]=febStruct2febFile(FEBioInputStruct)
% ------------------------------------------------------------------------
%
% This function uses the input FEBioInputStruct to generate an FEBio .feb file.
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
% 2018/01/17: Added support for febiospec 2.5 and made it default
% 2018/01/18: Decided to make 2.5 available in new seperate system
%------------------------------------------------------------------------

%%

warning('febStruct2febFile is depricated and will be removed in future releases. Update your codes to use febioStruct2xml, see HELP_febioStruct2xml (furthermore febio_spec version 2.5 is recommended)');

%%

dispStartTitleGibbonCode('Writing FEBio XML object');

%% Set default display setting if missing
if ~isfield(FEBioInputStruct,'disp_opt')
    FEBioInputStruct.disp_opt=0;
end

%% Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field

febioSpecVersion='2.0';

%Set febio_spec
febio_spec = domNode.getDocumentElement;
if ~isfield(FEBioInputStruct,'febio_spec')
    FEBioInputStruct.febio_spec.version=febioSpecVersion;
elseif ~isfield(FEBioInputStruct.febio_spec,'version')
    FEBioInputStruct.febio_spec.version=febioSpecVersion;
end
febio_spec.setAttribute('version',FEBioInputStruct.febio_spec.version); %Adding version attribute

disp(['Using febio_spec: ',FEBioInputStruct.febio_spec.version]);

%% Add comment if present
if isfield(FEBioInputStruct,'commentField')
    commentString = FEBioInputStruct.commentField;
else %Default comment
    commentString = ['Created using GIBBON, ',datestr(now)];
end
commentNode = domNode.createComment(commentString);
febio_spec.appendChild(commentNode);

%% DEFINING MODULE LEVEL
if ~isfield(FEBioInputStruct,'Module')
    FEBioInputStruct.Module.Type='solid'; %Use solid as default module
end
domNode=addModuleLevel_FEB(domNode,FEBioInputStruct);

%% DEFINE CONTROL SECTION
if isfield(FEBioInputStruct,'Control')
    domNode=addControlLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINING GLOBALS LEVEL
if isfield(FEBioInputStruct,'Globals')
    domNode=addGlobalsLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINING MATERIAL LEVEL
if isfield(FEBioInputStruct,'Materials')
    domNode=addMaterialLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINING GEOMETRY LEVEL

if isfield(FEBioInputStruct,'Geometry')
    writeMethod=1; %Use to switch between text mode (1) or XML model (2)
    switch writeMethod
        case 1 % TEXT FILE PARSING (faster for large models)
            domNode=addGeometryLevel_TXT(domNode,FEBioInputStruct);
        case 2 %XML PARSING
            domNode=addGeometryLevel_FEB(domNode,FEBioInputStruct);
    end
end

%% DEFINE BOUNDARY CONDITIONS LEVEL AND LOADS LEVEL
if isfield(FEBioInputStruct,'Boundary')    
    [domNode]=addBoundaryLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINE LOAD LEVEL
if isfield(FEBioInputStruct,'Loads')    
    [domNode]=addLoadsLevel_FEB(domNode,FEBioInputStruct);    
end

%% DEFINE CONTACT
if isfield(FEBioInputStruct,'Contact')
    [domNode]=addContactLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINE CONSTRAINTS LEVEL
if isfield(FEBioInputStruct,'Constraints')
    [domNode]=addConstraintsLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINE LOADDATA LEVEL
if isfield(FEBioInputStruct,'LoadData')
    [domNode]=addLoadDataLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINE STEP LEVEL
if isfield(FEBioInputStruct,'Step')
    [domNode]=addStepLevel_FEB(domNode,FEBioInputStruct);
end

%% DEFINE OUTPUT LEVEL
domNode=addOutputLevel_FEB(domNode,FEBioInputStruct);

%% CREATE OUTPUT OR EXPORT XML FILE

switch nargout
    case 0
        disp('Writing .feb file');
        xmlwrite_xerces(FEBioInputStruct.run_filename,domNode); %Custom XML write function 
%         if isfield(FEBioInputStruct,'topCommentLine')
%             exportFEB_XML(FEBioInputStruct.run_filename,domNode,FEBioInputStruct.topCommentLine); % Saving XML file
%         else
%             exportFEB_XML(FEBioInputStruct.run_filename,domNode); % Saving XML file
%         end
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
