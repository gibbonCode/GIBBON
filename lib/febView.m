function [varargout]=febView(varargin)

%% Parse input

switch nargin
    case 1
        xmlSpec=varargin{1};
        viewerOpt=1; 
    case 2
        xmlSpec=varargin{1};
        viewerOpt=varargin{2};
end

if ischar(xmlSpec)
    exportXML=0;
    fileName=xmlSpec;
else
    pathName=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','temp');
    fileName=fullfile(pathName,'temp.xml');
    if isstruct(xmlSpec) %Assuming Febio structure        
        febioStruct2xml(xmlSpec,fileName);
    else %Assuming xmlSpec is a domNode
        xmlwrite_xerces(fileName,xmlSpec); %Custom XML write function
    end
end

[hFig]=xmlView(fileName,viewerOpt);
varargout{1}=hFig;

end

