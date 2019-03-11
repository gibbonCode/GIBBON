function FEBioPath=getFEBioPath

% function FEBioPath=getFEBioPath
% ------------------------------------------------------------------------
% This function reads a user specified FEBio path name from a config text
% file. The text file may contain multiple paths (e.g. for multiple
% operational systems), however the first valid patch is always used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/10/05 Updated for GIBBON
% 2015/10/05 Enabled searching valid paths in multi-line config. file.
% 2017/06/10 changed empty output to an empty character array. 
%------------------------------------------------------------------------

%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths
try 
    [T]=txtfile2cell(fileName);    
catch
    T={};
end

%Get first valid path
FEBioPath='';
if ~isempty(T)
    for q=1:1:numel(T)
        pathName=T{q};
        if exist(pathName,'file')==2
            FEBioPath=pathName;
            return
        end
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
