clear; close all; clc;

docPath=fileparts(fileparts(mfilename('fullpath')));
htmlFolderName=fullfile(docPath,'html');
markDownFileName=fullfile(docPath,'Documentation.md');

%% Getting DICOM file names and image dimension parameters
files = dir(fullfile(htmlFolderName,'*.html'));
files={files(1:end).name};
files=sort(files(:));
NumberOfFiles=numel(files);

%%

T={...
    '---';...
    'layout: page';...
    'title: "Documentation"';...
    'logo: "img/home-bg.jpg"';...
    'description: "Help and demo links"';...
    'header-img: "img/home-bg.jpg"';...
    '---';...
    ' ';...    
    '__The documentation is a work in progress__. Not all functions have associated help files and not all functionality is covered by the demos. The help files currently cover about 50% of the functions and the demos mainly cover the use of FEBio.  ';...
    ' ';...    
    '#### Table of content';...
    '* [The MATLAB integrated help](#helpMatlab)';...          
    '* [Function help files](#help)';...
    '* [Demo files](#demo)';...          
    ' ';...
    '## The MATLAB integrated help <a name="helpMatlab"></a>  ';...
    'Follow the installation instructions to integrate and access GIBBON documentation from within MATLAB. The name for all function help files (the files that generate the help/documentation when published with MATLAB) starts with `HELP_`, all demo files start with `DEMO_`. This way users may explore/open/edit these files by typing `open HELP_functionName` or `open DEMO_functionName` in the MATLAB command window';...
    ' ';...
    };

T{end+1}=['## Function help files <a name="help"></a>   ' ];
T{end+1}='';
for q=1:1:NumberOfFiles
    [~,fileNameNow,~]=fileparts(files{q});
    if strfind(files{q},'HELP_')        
        fileNamePNG=fullfile(htmlFolderName,[fileNameNow,'.png']);
        if exist(fileNamePNG,'file')            
           T{end+1}=['[',fileNameNow,'](html/',files{q},') ','![',fileNameNow,'](html/',[fileNameNow,'.png'],'){:height="40px"}  ']; 
        else
            T{end+1}=['[',fileNameNow,'](html/',files{q},')  '];        
        end        
    end
end

T{end+1}=['## Demo files <a name="demo"></a>  '];
T{end+1}='';
for q=1:1:NumberOfFiles
    [~,fileNameNow,~]=fileparts(files{q});
    if strfind(files{q},'DEMO_')
        if exist(fileNamePNG,'file')            
           T{end+1}=['[',fileNameNow,'](html/',files{q},') ','![',fileNameNow,'](html/',[fileNameNow,'.png'],'){:height="40px"}  ']; 
        else
            T{end+1}=['[',fileNameNow,'](html/',files{q},')  '];        
        end        
    end
end

cell2txtfile(markDownFileName,T,0);
 
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
