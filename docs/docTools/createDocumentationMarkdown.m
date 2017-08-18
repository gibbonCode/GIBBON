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
        if exist(fileNamePNG,'file')==2
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
        fileNamePNG=fullfile(htmlFolderName,[fileNameNow,'.png']);
        if exist(fileNamePNG,'file')==2
            T{end+1}=['[',fileNameNow,'](html/',files{q},') ','![',fileNameNow,'](html/',[fileNameNow,'.png'],'){:height="40px"}  '];
        else
            T{end+1}=['[',fileNameNow,'](html/',files{q},')  '];
        end
    end
end

cell2txtfile(markDownFileName,T,0);

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
