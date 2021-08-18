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
