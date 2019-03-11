function write_XML_no_extra_lines(fileName,domNode)

% function write_XML_no_extra_lines(fileName,domNode)
% ------------------------------------------------------------------------
%
% This function writes XML files
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log:
% 2014/05/27: Updated for GIBBON
% 2017/06/26: Implemented xmlwrite_xerces to avoid bug in MATLAB 2017. 
%------------------------------------------------------------------------

%%

%Write to text file
%xmlwrite(fileName,domNode); %MATLAB's xmlwrite extremely slow for v2017a-2017b

xmlwrite_xerces(fileName,domNode); %Custom XML write function

%Import back into cell array
[T]=txtfile2cell(fileName);

%Save to txt file while skipping "empty lines"
cell2txtfile(fileName,T,1);
 
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
