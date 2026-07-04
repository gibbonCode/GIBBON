function [testFileList]=getTestFiles(testSet)

% function [testFileList]=getTestFiles(testSet)
% ------------------------------------------------------------------------
% This function obtains the list of test files from the documentation
% folder. The test files are either the HELP_* files (testSet='help'), the
% DEMO_* files (testSet='demo') or all (testSet='all'). 
%
% Change log:
% 2021/08/05 Created
% ------------------------------------------------------------------------

%%
testFileList=fullfile(fileparts(fileparts(mfilename('fullpath'))),'docs');

%Get list of M-files in documentation folder
allFiles_publish = dir(fullfile(testFileList,'*.m'));
allFiles_publish={allFiles_publish(1:end).name};
allFiles_publish=sort(allFiles_publish(:));

logicExclude=~gcontains(allFiles_publish,'obj_DEMO_'); %Skip objective functions
logicDemo=gcontains(allFiles_publish,'DEMO_'); %Logic for demo files
logicHelp=gcontains(allFiles_publish,'HELP_'); %Logic for help files

%Get list of files, either all or just help or demo files
switch lower(testSet) %Lower so both HELP and help etc work
    case 'all'
        testFileList=allFiles_publish((logicHelp | logicDemo) & logicExclude);
    case 'demo'
        testFileList=allFiles_publish(logicDemo & logicExclude);
    case 'help'
        testFileList=allFiles_publish(logicHelp & logicExclude);
    otherwise
        error([testMode,' is not a valid option. Use: all, demo, or help']);
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
