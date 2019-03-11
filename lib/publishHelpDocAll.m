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
