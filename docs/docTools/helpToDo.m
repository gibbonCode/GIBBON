clear; close all; clc; 

%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(fileparts(filePath)));

libFolder=fullfile(toolboxPath,'lib');
helpFolder=fullfile(toolboxPath,'docs');
outputFolder=fullfile(toolboxPath,'docs','todo');

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

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
