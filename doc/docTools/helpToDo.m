clear; close all; clc; 

%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(fileparts(filePath)));

libFolder=fullfile(toolboxPath,'lib');
helpFolder=fullfile(toolboxPath,'doc');
outputFolder=fullfile(toolboxPath,'doc','todo');

%%

%Get list of M-files in library
fileListLib = dir(fullfile(libFolder,'*.m'));
fileListLib={fileListLib(1:end).name};
fileListLib=sort(fileListLib(:));
NumberOfLibFiles=numel(fileListLib);

%Get list of M-files in help folder
fileListHelp = dir(fullfile(helpFolder,'*.m'));
fileListHelp={fileListHelp(1:end).name};
fileListHelp=sort(fileListHelp(:));
NumberOfHelpFiles=numel(fileListHelp);

%%

if exist(outputFolder,'dir')==0
   mkdir(outputFolder);  
end
addpath(outputFolder);

for q=1:1:NumberOfLibFiles
   libFileName=fileListLib{q};
   libFileNameFull=fullfile(libFolder,libFileName);
   helpFileName=fullfile(helpFolder,['HELP_',libFileName]);
   helpFileNameNew=fullfile(outputFolder,['HELP_',libFileName]);
   
   
   if exist(helpFileName,'file')==0
       funcName=libFileName(1:end-2);
       [T_func]=txtfile2cell(libFileNameFull); %A cell containing the M-file text
       funcLine=T_func{1}; %First line containing function statement
       funcLine = strtrim(funcLine); %First line with leading/trailing spaces removed
       funcLine = regexprep(funcLine,' +',' '); %Remove possible extra spaces
       codeSyntaxSimple=funcLine(10:end);
              
       %Build M-file text
       T={};
       T{1,1}=['%% ',funcName];
       T{end+1,1}=['% Below is a demonstration of the features of the |',funcName,'| function'];
       T{end+1,1}='';
       T{end+1,1}='%%';
       T{end+1,1}='clear; close all; clc;';
       T{end+1,1}='';
       T{end+1,1}='%% Syntax';
       T{end+1,1}=['% |',codeSyntaxSimple,';|'];
       T{end+1,1}='';
       T{end+1,1}='%% Description ';
       T{end+1,1}='% UNDOCUMENTED ';
       T{end+1,1}='%% Examples ';
       T{end+1,1}='% ';
       T{end+1,1}='%%';
       T{end+1,1}='% ';
       T{end+1,1}='% <<gibbVerySmall.gif>>';
       T{end+1,1}='% ';
       T{end+1,1}='% _*GIBBON*_ ';
       T{end+1,1}='% <www.gibboncode.org>';
       T{end+1,1}='% ';
       T{end+1,1}='% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>';
       
       cell2txtfile(helpFileNameNew,T,0);
       
%        publish(helpFileNameNew);

   end
end

