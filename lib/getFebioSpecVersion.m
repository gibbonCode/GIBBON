function [febio_spec]=getFebioSpecVersion(febXML)

% function [febio_spec]=getFebioSpecVersion(febXML)
% ------------------------------------------------------------------------
% Reads in the version attribute of the febio_spec field
% E.g. the XML entry:
% <febio_spec version="2.0">
% Would result in the output: febio_spec='2.0'
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/10/10
%------------------------------------------------------------------------

if ischar(febXML); %If input is a string assume is the filename for the XML
   febXML=xmlread(febXML);  
end

febio_spec=febXML.getElementsByTagName('febio_spec').item(0).getAttribute('version').toCharArray()';