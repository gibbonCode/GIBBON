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
 
%% <-- GIBBON footer text --> 
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
