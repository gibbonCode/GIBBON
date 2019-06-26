function changeFileExtensions(pathName,extOld,extNew)

% function changeFileExtensions(pathName,extOld,extNew)
%------------------------------------------------------------------------
% 
%
%
%------------------------------------------------------------------------

%%
if isempty(extOld)
    files = dir(pathName);
    files={files(1:end).name};
else
    files = dir(fullfile(pathName,['*.',extOld]));
    files={files(1:end).name};
end

for q=1:1:numel(files)
    oldNameFull=fullfile(pathName,files{q});
    if exist(oldNameFull,'file')
        [~,oldName,extFile] = fileparts(oldNameFull);
        if isempty(extNew)
            newNameFull=fullfile(pathName,oldName);
        else
            newNameFull=fullfile(pathName,[oldName,'.',extNew]);
        end
        if strcmp(oldNameFull,newNameFull)==0
            if isempty(extOld)
                if strcmp(extFile,extOld)
                    movefileNow(oldNameFull,newNameFull);
                end
            else                
                movefileNow(oldNameFull,newNameFull);
            end
        end
    end
end

end

function movefileNow(oldNameFull,newNameFull)
    movefile(oldNameFull,newNameFull);
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
