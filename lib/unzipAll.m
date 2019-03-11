function unzipAll(varargin)

% function unzipAll(pathName,pathNameOutput)
% ------------------------------------------------------------------------
%
% This function unzips all zip files in the folder pathName. The first
% input is the folder path name, the following inputs are possible inputs
% of the unzip function. 
%
% See also: dicomanon
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/09: Created
%------------------------------------------------------------------------

%%
switch nargin 
    case 1
        pathName=varargin{1};
        varOutput=pathName;
    case 2
        pathName=varargin{1};
        varOutput=varargin{2};
end
%%

if ischar(varOutput) %Process only input folder
    files = dir(fullfile(pathName,'*.zip'));
    files={files(1:end).name};
    files=sort(files(:));
    if ~isempty(files)
        if numel(files)>1
            hw = waitbar(0,'Unzipping files...');
            for q=1:1:numel(files)
                fileName=fullfile(pathName,files{q});
                unzip(fileName,varOutput);
                waitbar(q/numel(files),hw,['Unzipping files...',num2str(round(100.*q/numel(files))),'%']);
            end
            close(hw)
        else
            fileName=fullfile(pathName,files{1});
            unzip(fileName,varOutput);
        end
    end
else %Process input folder and all subfolders
    %Process current
    unzipAll(pathName,pathName); 
    
    %Process sub-folders
    [pathNames]=getSubPaths(pathName);
    if ~isempty(pathNames)
        for q=1:1:numel(pathNames)
            pathNameNow=pathNames{q};
            unzipAll(pathNameNow,pathNameNow);
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
