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
if nargout>0
    varargout{1}=hFig;
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
