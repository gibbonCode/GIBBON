function [varargout]=febStruct2febFile_v1p2(FEB_struct)

% function [varargout]=febStruct2febFile_v1p2(FEB_struct)
% ------------------------------------------------------------------------
%
% This function uses the input FEB_struct to generate a basic .feb
% file suitable for import into PreView. 
%
% 
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/05/27
%------------------------------------------------------------------------

%%

dispStartTitleGibbonCode('Writing FEBio XML object');

%% Initialize docNode object
docNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field

%Set febio_spec
febio_spec = docNode.getDocumentElement;
if ~isfield(FEB_struct,'febio_spec')
    FEB_struct.febio_spec.version='1.2';
elseif ~isfield(FEB_struct.febio_spec,'version')
    FEB_struct.febio_spec.version='1.2';
end
febio_spec.setAttribute('version',FEB_struct.febio_spec.version); %Adding version attribute 

%% DEFINING MODULE LEVEL
docNode=addModuleLevel_FEB_v1p2(docNode,FEB_struct);

%% DEFINE CONTROL SECTION
if isfield(FEB_struct,'Control');
    docNode=addControlLevel_FEB_v1p2(docNode,FEB_struct);    
end

%% DEFINING GLOBALS LEVEL
docNode=addGlobalsLevel_FEB_v1p2(docNode,FEB_struct);

%% DEFINING MATERIAL LEVEL
docNode=addMaterialLevel_FEB_v1p2(docNode,FEB_struct);

%% DEFINING GEOMETRY LEVEL
writeMethod=1;
switch writeMethod
    case 1 % TEXT FILE PARSING (faster for large arrays)
                docNode=addGeometryLevel_TXT_v1p2(docNode,FEB_struct);
    case 2 %XML PARSING
                docNode=addGeometryLevel_FEB_v1p2(docNode,FEB_struct);
end

%% DEFINE OUTPUT LEVEL
docNode=addOutputLevel_FEB_v1p2(docNode,FEB_struct);

%% DEFINE BOUNDARY CONDITIONS LEVEL AND LOADS LEVEL
if isfield(FEB_struct,'Boundary')
    [docNode]=addBoundaryLevel_FEB_v1p2(docNode,FEB_struct);
end

%% DEFINE LOAD LEVEL
if isfield(FEB_struct,'Loads')
    [docNode]=addLoadsLevel_FEB_v1p2(docNode,FEB_struct);
end

%% DEFINE LOADDATA LEVEL
if isfield(FEB_struct,'LoadData')
    [docNode]=addLoadDataLevel_FEB_v1p2(docNode,FEB_struct);
end

%% DEFINE CONSTRAINTS LEVEL
if isfield(FEB_struct,'Constraints');
    [docNode]=addConstraintsLevel_FEB_v1p2(docNode,FEB_struct);
end

%% DEFINE STEP LEVEL
if isfield(FEB_struct,'Step');
    [docNode]=addStepLevel_FEB_v1p2(docNode,FEB_struct);
end

%% CREATE OUTPUT OR EXPORT XML FILE 

switch nargout
    case 0
        disp('Writing .feb file');        
        if isfield(FEB_struct,'topCommentLine')            
            exportFEB_XML(FEB_struct.run_filename,docNode,FEB_struct.topCommentLine); % Saving XML file
        else
            exportFEB_XML(FEB_struct.run_filename,docNode); % Saving XML file
        end                
    case 1
        varargout{1}=docNode;
end
dispDoneGibbonCode;

end


