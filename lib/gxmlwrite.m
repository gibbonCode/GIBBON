function [varargout]=gxmlwrite(varargin)

% function [varargout]=gxmlwrite(fileName,domNode,writeMode)
% ------------------------------------------------------------------------
% This function serializes an XML Document Object Model node (domNode) to a
% file (fileName). The XML file is created using either the Xerces parser
% (writeMode=1) or MATLAB's parser (writeMode=2). The MATLAB parses is
% known to have issues for MATLAB 2017a-b. 
%
% The function is similar to MATLAB's xmlwrite function. However by default
% it uses the Xerces parser rather than MATLAB's parser. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log: 
% 2017/08/18: Created
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1        
        domNode=varargin{1};
        fileName=[];
        writeMode=2;
    case 2
        fileName=varargin{1};
        domNode=varargin{2};
        writeMode=1;
    case 3
        fileName=varargin{1};
        domNode=varargin{2};
        writeMode=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

%% Write XML file if needed

if ~isempty(fileName)
    switch writeMode
        case 1 %Xerces parser
            xmlwrite_xerces(fileName,domNode); 
        case 2 %MATLAB's parse
            xmlwrite(fileName,domNode); 
    end
end

%% Formulate output if needed
if nargout>0
    %Generate string
    switch writeMode
        case 1 %Xerces parser
            XML_str = xmlwrite_xerces(domNode); 
        case 2 %MATLAB's parse
            XML_str = xmlwrite(domNode); 
    end
    varargout{1}=XML_str;
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
