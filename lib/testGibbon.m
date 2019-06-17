function testGibbon(varargin)

% function testGibbon(testSet,testMode,approveQuestion,startInd)
%-------------------------------------------------------------------------
% This function can be used to run, test, and if desired publish, all HELP_
% and DEMO_ files. These files also define the help and documentation and
% therefore publishing these codes generates the integrated help files.
%
%   testGibbon('all','test',0,1); %Test run all files   
%   testGibbon('help','test',0,1); %Test run all HELP_ files
%   testGibbon('demo','test',0,1); %Test run all DEMO_ files
%   testGibbon('all','pub',0,1); %Test run and publish all files
%   testGibbon('help','pub',0,1); %Test run and publish all HELP_ files
%   testGibbon('demo','pub',0,1); %Test run and publish all DEMO_ files
%-------------------------------------------------------------------------

%% Parse input
switch nargin
    case 0
        testSet='all';
        testMode='test';
        approveQuestion=0;
        startLoc=1;
    case 1
        testSet=varargin{1};
        testMode='test';
        approveQuestion=0; 
        startLoc=1;
    case 2
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=0;
        startLoc=1;
    case 3
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=varargin{3};    
        startLoc=1;
    case 4
        testSet=varargin{1};
        testMode=varargin{2};
        approveQuestion=varargin{3};
        startLoc=varargin{4};
end

%% Running tests, publishing codes

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
    otherwise
        error([testMode,' is not a valid option. Use: all, demo, or help']);
end

%% Get start index
if ischar(startLoc) || isstring(startLoc)
    startInd=find(gcontains(testFileList,startLoc));
    if isempty(startInd)
        error(['Start file "',startLoc,'" not found in list'])
    end
else
    startInd=startLoc;
end

%%

addpath(testFolder);
cd(testFolder); %Make testFolder current directory

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
            try
                run(mFileNow);
                drawnow;
                pause(1);
            catch ME
                load('tempPub.mat'); %Load variables
                disp(['Error in: ',mFileNow]);
                disp(['To restart from this file use start index: ',num2str(q_test)]);
                disp(' ');
                rethrow(ME);
            end            
        case 'pub' %Publish
            try
                gpublish(mFileNow);
                % publish(mFileNow,'catchError',false,'figureSnapMethod','getframe','maxWidth',1000);
                drawnow;
                pause(1);
            catch ME                             
                load('tempPub.mat'); %Load variables
                disp(['Error in: ',mFileNow]);
                disp(['To restart from this file use start index: ',num2str(q_test)]);
                disp(' ');
                rethrow(ME);
            end
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
