function testGibbon(varargin)

switch nargin
    case 0
        testSet='all';
        testMode='test';
        approveQuestion=1;
        startInd=1;
    case 1
        testSet=varargin{1};
        testMode='test';
        approveQuestion=1; 
        startInd=1;
    case 2
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=1;
        startInd=1;
    case 3
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=varargin{3};    
        startInd=1;
    case 4
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=varargin{3};
        startInd=varargin{4};
end

%%

testFolder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'docs');

%Get list of M-files in documentation folder
allFiles_publish = dir(fullfile(testFolder,'*.m'));
allFiles_publish={allFiles_publish(1:end).name};
allFiles_publish=sort(allFiles_publish(:));
logicNot=~cellfun(@isempty,strfind(allFiles_publish,'obj_DEMO_')); %Skip objective functions

logicDemo=~cellfun(@isempty,strfind(allFiles_publish,'DEMO_')); %Logic for demo files
logicHelp=~cellfun(@isempty,strfind(allFiles_publish,'HELP_')); %Logic for help files

%Get list of files, either all or just help or demo files
switch testSet
    case 'all'
        testFileList=allFiles_publish((logicHelp | logicDemo) & ~logicNot);
    case 'demo'
        testFileList=allFiles_publish(logicDemo & ~logicNot);
    case 'help'
        testFileList=allFiles_publish(logicHelp & ~logicNot);        
end

%%

addpath(testFolder);

for q_test=startInd:1:numel(testFileList)

    fileMessage=['testGibbon -> Test file:',num2str(q_test),' of ',num2str(numel(testFileList)),' ',testFileList{q_test}];
    disp(' ');
    disp(fileMessage);
    disp(' ');
    
    mFileNow=fullfile(testFolder,testFileList{q_test});
    initialVars_publish = who;        
    save('tempPub.mat'); %Save current variables
    
    switch testMode
        case 'test' %Test
            run(mFileNow);
        case 'pub' %Publish
            publish(mFileNow,'catchError',false,'figureSnapMethod','getframe','maxHeight',800);
    end
    
    load('tempPub.mat'); %Load variables
    delete('tempPub.mat'); %Clean up
   
    if approveQuestion
        choice = questdlg([fileMessage,'. Done. Do you want to proceed?'],testFileList{q_test},'Yes','No','Yes');        
        switch choice
            case 'Yes'
                
            case 'No'
                edit(mFileNow);
                break
        end
    end
    clearvars('-except',initialVars_publish{:});
    close all;
    
end
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
