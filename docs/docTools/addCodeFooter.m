clear; close all; clc;

%%

toolboxPath=fileparts(fileparts(fileparts(mfilename('fullpath'))));
docPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs');
libPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'lib');
docToolsPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs','docTools');

pathNames={docToolsPath,libPath,docPath};

%% Prepare boilerplate

licenseBoilerPlate=fullfile(toolboxPath,'licenseBoilerPlate.txt');
[T]=txtfile2cell(licenseBoilerPlate);
size(T)

footerTargetText='%% <-- GIBBON footer text --> ';

%Add comment symbols
for q=1:1:numel(T)
    T{q}=['% ',T{q}];
end

%Add target header
T=[{footerTargetText};T(1:end)];

%%

fileExtension='.m';

%%

numPaths=numel(pathNames); 
for q_path=1:1:numPaths
    pathName=pathNames{q_path};     
    files = dir(fullfile(pathName,'*.m'));
    files={files(1:end).name};
    files=sort(files(:));
    numFiles=numel(files);    
    for q_file=1:1:numFiles
        fileName=fullfile(pathName,files{q_file});
        [T_now]=txtfile2cell(fileName);
        
        T_now(end+1:end+numel(T))=T;
        cell2txtfile(fileName,T_now,0);
        
    end
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
