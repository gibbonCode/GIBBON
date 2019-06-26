function [varargout]=ckmkdir(dirName)

% function [varargout]=ckmkdir(dirName)
% ------------------------------------------------------------------------
% Check for the existance of the directory dirName and make the directory
% if it does not. The optional output is the exitFlag which is 0 if it did
% not exist and 1 if it did. 
%
% ------------------------------------------------------------------------

%%

% Check for existence of directory
existFlag=exist(dirName,'dir');
if ~existFlag %Check if folder exists   
    mkdir(dirName); %Make the directory if it does not exist
end

%Store output
if nargout==1
    varargout{1}=existFlag;
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
