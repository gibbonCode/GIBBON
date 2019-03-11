function setFEBioPath(FEBioPathSpec)

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
            end
        end
    otherwise
        if exist(FEBioPathSpec,'file')==2
            T{2}=FEBioPathSpec;
        end
end

cell2txtfile(fileName,T,0,0); %Write cell to config file

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
