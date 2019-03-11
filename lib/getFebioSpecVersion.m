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
