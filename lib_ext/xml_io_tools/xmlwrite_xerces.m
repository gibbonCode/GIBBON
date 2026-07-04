function varargout=xmlwrite_xerces(varargin)

% ----------------------------------------------------------
% function varargout=xmlwrite_xerces(FILENAME,DOMNODE)
% XMLWRITE_XERCES Serialize an XML Document Object Model node using Xerces parser. 
%  xmlwrite_xerces(FILENAME,DOMNODE) serializes the DOMNODE to file FILENAME.
%
% The function xmlwrite_xerces is very similar the Matlab function xmlwrite 
% but works directly with the XERCES java classes (written by Apache XML 
% Project) instead of the XMLUtils class created by Mathworks. Xerces files
% are provided in standard MATLAB instalation and live in root\java\jarext
% directory. 
%
% This function is needed as a work-around for a bug in XMLUtils library
% which can not write CDATA SECTION nodes correctly. Also Xerces and 
% XMLUtils libraries handle namespaces differently.  
%
% Examples:
%   % See xmlwrite examples this function have almost identical behavior.
%  
% Advanced use:
%  FILENAME can also be a URN, java.io.OutputStream or java.io.Writer object
%  SOURCE can also be a SAX InputSource, JAXP Source, InputStream, or 
%    Reader object
%
% Change log:
% 2017/06/26 This code was obtained from 'xml_io_tools' copyright (c) 2007,
% Jaroslaw Tuszynski (BSD-2-Clause license) 
% (https://nl.mathworks.com/matlabcentral/fileexchange/12907-xml-io-tools). 
% Although the original code is labelled as having been written by Alberto Amaro
% (2007/02/22). 
% 2017/06/26 Kevin Moerman fixed code to work with MATLAB 2016,2017. This
% involved changing a setting from 'on' (string) to 1 (scalar). 
% 2017/06/26 Kevin Moerman determined that this code fixes bug in 2017's
% xmlwrite. i.e.
% javaMethod('serializeXML','com.mathworks.xml.XMLUtils',source,result); is
% very slow in MATLAB 2017a-b.
%
% To do: 
% 2017/08/18 Fix error that occurs for string output with singular input
% ----------------------------------------------------------

%%

returnString = false;
if length(varargin)==1
    returnString = true;
    result = java.io.StringWriter;
    source = varargin{1};
else
    result = varargin{1};
    if ischar(result)
      % Using the XERCES classes directly, is not needed to modify the
      % filename string. So I have commented this next line
      %  result = F_xmlstringinput(result,false);
    end
    
    source = varargin{2};
    if ischar(source)
        source = F_xmlstringinput(source,true);
    end
end

% SERIALIZATION OF THE DOM DOCUMENT USING XERCES CLASSES DIRECTLY

% 1) create the output format according to the document definitions
% and type
objOutputFormat = org.apache.xml.serialize.OutputFormat(source);
set(objOutputFormat,'Indenting',1); %KMM: Changed from set(objOutputFormat,'Indenting','on');

% 2) create the output stream. In this case: an XML file
objFile = java.io.File(result);
objOutputStream = java.io.FileOutputStream(objFile);

% 3) Create the Xerces Serializer object
objSerializer= org.apache.xml.serialize.XMLSerializer(objOutputStream,objOutputFormat);

% 4) Serialize to the XML files
javaMethod('serialize',objSerializer,source);

% 5) IMPORTANT! Delete the objects to liberate the XML file created
objOutputStream.close;

if returnString
    varargout{1}=char(result.toString);
end

%% ========================================================================
 function out = F_xmlstringinput(xString,isFullSearch,varargin)
% The function F_xmlstringinput is a copy of the private function:
% 'xmlstringinput' that the original xmlwrite function uses.

if isempty(xString)
    error('Filename is empty');
elseif ~isempty(findstr(xString,'://'))
    %xString is already a URL, most likely prefaced by file:// or http://
    out = xString;
    return;
end

xPath=fileparts(xString);
if isempty(xPath)
    if nargin<2 || isFullSearch
        out = which(xString);
        if isempty(out)
            error('xml:FileNotFound','File %s not found',xString);
        end
    else
        out = fullfile(pwd,xString);
    end
else
    out = xString;
    if (nargin<2 || isFullSearch) && ~exist(xString,'file')
        %search to see if xString exists when isFullSearch
        error('xml:FileNotFound','File %s not found',xString);
    end
end

%Return as a URN
if strncmp(out,'\\',2)
    % SAXON UNC filepaths need to look like file:///\\\server-name\
    out = ['file:///\',out];
elseif strncmp(out,'/',1)
    % SAXON UNIX filepaths need to look like file:///root/dir/dir
    out = ['file://',out];
else
    % DOS filepaths need to look like file:///d:/foo/bar
    out = ['file:///',strrep(out,'\','/')];
end

