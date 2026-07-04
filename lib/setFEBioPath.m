function setFEBioPath(FEBioPathSpec)

% function setFEBioPath(FEBioPathSpec)
% -----------------------------------------------------------------------
% This function sets the FEBio path defined by FEBioPathSpec. Setting the
% path means it is stored in the GIBBON configuration file:
% /GIBBON/config/FEBioPath.txt. 
% 
% -----------------------------------------------------------------------

%%
filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths

T{1}='#List paths for febio here, multiple can be given e.g. for use on several platforms. The first valid path found is always used.';

%Set paths in cell
switch class(FEBioPathSpec)
    case 'cell'
        for q=1:1:numel(FEBioPathSpec)
            pathNameNow=FEBioPathSpec{q};
            if exist(pathNameNow,'file')==2
                T{q+1}=FEBioPathSpec{q};
            else
                error('Path specified does not exist')
            end
        end
    otherwise
        if exist(FEBioPathSpec,'file')==2
            T{2}=FEBioPathSpec;
        else
            error('Path specified does not exist')
        end
end

cell2txtfile(fileName,T,0,0); %Write cell to config file

if contains(lower(FEBioPathSpec),'febio2')
        warning('FEBio2 detected. FEBio2 support is depricated. Please upgrade to FEBio3'); 
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
