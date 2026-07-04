function [pathNames]=getSubPaths(pathName, ignoreHidden)

% function [pathNames]=getSubPaths(pathName, ignoreHidden)
% ------------------------------------------------------------------------
%
% This function (based on the GENPATH command) creates the output pathNames
% which contains all folders and sub-folders within the folder specified by
% the input pathName.
%
% If the optional input ignoreHidden is true, folders starting with a dot
% will be ignored.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2013/04/18 Created
% 2017/06/01 Fixed bug in relation to operational system differences
% 2026/06/05 Added option to ignore hidden folders
%------------------------------------------------------------------------

%%
if ispc % semi-colon splits paths
    strPattern=[filesep,';'];    
else %colon splits paths for Linux and mac
    strPattern=':';        
end
pathNames=regexp(genpath(pathName),strPattern, 'split');
pathNames=pathNames(2:end-1)';

if nargin > 1 && ignoreHidden
    pathNames(cellfun(@isHiddenFolder,pathNames)) = [];
end

end

function yn = isHiddenFolder(path)
    assert(isfolder(path))
    [~, name] = fileparts(path);
    yn = isempty(name) || startsWith(name, '.');
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
