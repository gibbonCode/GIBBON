function varargout=xmlwriteCustom(varargin)
%XMLWRITE  Serialize an XML Document Object Model node.
%   XMLWRITE(FILENAME,DOMNODE) serializes the DOMNODE to file FILENAME.
%
%   S = XMLWRITE(DOMNODE) returns the node tree as a string.
%
%   Example:
%   % Create a sample XML document.
%   docNode = com.mathworks.xml.XMLUtils.createDocument('root_element')
%   docRootNode = docNode.getDocumentElement;
%   docRootNode.setAttribute('attribute','attribute_value');
%   for i=1:20
%      thisElement = docNode.createElement('child_node');
%      thisElement.appendChild(docNode.createTextNode(sprintf('%i',i)));
%      docRootNode.appendChild(thisElement);
%   end
%   docNode.appendChild(docNode.createComment('this is a comment'));
%
%   % Save the sample XML document.
%   xmlFileName = [tempname,'.xml'];
%   xmlwrite(xmlFileName,docNode);
%   edit(xmlFileName);
%
%   See also XMLREAD, XSLT.

%   Copyright 1984-2006 The MathWorks, Inc.

%    Advanced use:
%       FILENAME can also be a URN, java.io.OutputStream or
%                java.io.Writer object
%       SOURCE can also be a SAX InputSource, JAXP Source,
%              InputStream, or Reader object

% This is the XML that the help example creates:
% <?xml version="1.0" encoding="UTF-8"?>
% <root_element>
%     <child_node>1</child_node>
%     <child_node>2</child_node>
%     <child_node>3</child_node>
%     <child_node>4</child_node>
%     ...
%     <child_node>18</child_node>
%     <child_node>19</child_node>
%     <child_node>20</child_node>
% </root_element>
% <!--this is a comment-->

filename = [];

returnString = false;
if length(varargin)==1
    returnString = true;
    result = java.io.StringWriter;
    source = varargin{1};
else
    result = varargin{1};
    if ischar(result)
        filename = result;
        result = xmlstringinput(result,false);
        % This strips off the extra stuff in the resolved file.  Then,
        % we are going to use java to put it in the right form.
        if strncmp(result, 'file:', 5)
           result = regexprep(result, '^file:///(([a-zA-Z]:)|[\\/])','$1');
           result = strrep(result, 'file://', '');
           temp = java.io.File(result);
           result = char(temp.toURI());
        end
    elseif ~isa(result, 'java.io.Writer') && ~isa(result, 'java.io.OutputStream')
            error(message('MATLAB:xmlwrite:IncorrectFilenameType'));
    end
    
    source = varargin{2};
    if ischar(source)
        source = xmlstringinput(source,true);
    end
end

% The JAXP-approved way to serialize a 
% document is to run a null transform.
% This is a JAXP-compliant static convenience method
% which does exactly that.
javaMethod('serializeXML',...
    'com.mathworks.xml.XMLUtils',...
    source,result);

if returnString
%     varargout{1}=char(result.toString);
varargout{1}=(result.toString); %MODIFIED KMM 2014/01/14
else
    %this notifies the operating system of a file system change.  This
    %probably doesn't work if the user passed in the filename in the form
    %of file://filename, but it would probably be more trouble than it is
    %worth to resolve it.  It should be harmless in that case.
    if ischar(result) && strncmp(result, 'file:', 5)
        fschange(fileToDirectory(filename));
    end
end


function final_name = fileToDirectory(orig_name)
% This is adequate to resolve the full path since the call above to xmlstringinput 
% does not search the path when looking to write the file.
temp = fileparts(orig_name);
if isempty(temp)
    final_name = pwd;
else
    final_name = temp;
end