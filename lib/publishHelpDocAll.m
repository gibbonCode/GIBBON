function publishHelpDocAll(varargin)

switch nargin
    case 0 
        runHelp=1;
        runDemo=1;
    case 1
        runHelp=varargin{1};
        runDemo=1;
    case 2
        runHelp=varargin{1};
        runDemo=varargin{2};
end

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
helpPath=fullfile(toolboxPath,'docs');

%%

%Get list of M-files in help folder
allFiles_publish = dir(fullfile(helpPath,'*.m'));
allFiles_publish={allFiles_publish(1:end).name};
allFiles_publish=sort(allFiles_publish(:));

if runHelp
    logicHelp_publish=~cellfun(@isempty,strfind(allFiles_publish,'HELP_'));
else
    logicHelp_publish=false(size(allFiles_publish));
end

if runDemo
    logicDemo_publish=~cellfun(@isempty,strfind(allFiles_publish,'DEMO_'));
else
    logicDemo_publish=false(size(allFiles_publish));
end

logicNot=~cellfun(@isempty,strfind(allFiles_publish,'obj_DEMO_'));

pubFileList_publish=allFiles_publish((logicHelp_publish | logicDemo_publish) & ~logicNot);

numPubpFiles_publish=numel(pubFileList_publish);

%%

approveQuestion=1; 
catchErrorOpt=true;%false;
addpath(helpPath);

for q_publish=1:1:numPubpFiles_publish
    initialVars_publish = who;        
    disp(q_publish);
    save('tempPub.mat');
    pause(1);
    publish(fullfile(helpPath,pubFileList_publish{q_publish}),'catchError',catchErrorOpt,'figureSnapMethod','getframe','maxHeight',800);
    pause(1);
    load('tempPub.mat');
    disp(q_publish);
    
    disp(pubFileList_publish{q_publish});
   
    if approveQuestion
        choice = questdlg([pubFileList_publish{q_publish},' done. Proceed?'],pubFileList_publish{q_publish},'Yes','No','Yes');
        
        switch choice
            case 'Yes'
                
            case 'No'
                edit(fullfile(helpPath,pubFileList_publish{q_publish}));
                break
        end
    end
    clearvars('-except',initialVars_publish{:});
    disp(q_publish);
    close all;
end
%%
delete('tempPub.mat');

 
%% 
% ********** _license boilerplate_ **********
% 
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
