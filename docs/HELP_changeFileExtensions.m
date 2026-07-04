%% changeFileExtensions
% Below is a demonstration of the features of the |changeFileExtensions| function

%%
clear; close all; clc;

%% Syntax
% |changeFileExtensions(pathName,extOld,extNew);|

%% Description
% The |changeFileExtensions| function changes the file extensions for all
% files in the folder pathName, and with the extension extOld, to have the
% extension extNew.

%% Examples
%

%%
%Create example files in the data/temp directory. Here a set of files with
%the txt extension are created. Later these are changed to have a csv
%file extension.

%Create .txt files
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','temp','renameCheck');
mkdir(pathName)

n=5;
for q=1:1:n %Create n files
    fileID=fopen(fullfile(pathName,['temp',num2str(q),'.txt']),'w');
    fprintf(fileID,'%d\n',pi);
    fclose(fileID);
end

% Add one additional file with a different extension to show this file will
% not be effected. 
fileID=fopen(fullfile(pathName,'temp.txp'),'w');
fprintf(fileID,'%d\n',pi);
fclose(fileID);

%%
% Show current folder content
disp('Old folder content:')
ls(pathName)

%%
% Change file extensions
extOld='txt'; %Old extension
extNew='csv'; %New extension
changeFileExtensions(pathName,extOld,extNew)

%%
% Show current folder content
disp('New folder content:')
ls(pathName)

%%
% remove the temporary folder created for this example
rmdir(pathName,'s')

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
%%
% _*GIBBON footer text*_
%
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
%
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
%
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
