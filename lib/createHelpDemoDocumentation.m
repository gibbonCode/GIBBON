function createHelpDemoDocumentation

%% Get toolbox location

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
helpPath=fullfile(toolboxPath,'docs','html');

%% FIND HTML FILES

allFiles = dir(fullfile(helpPath,'*.html'));
allFiles={allFiles(1:end).name};
allFiles=sort(allFiles(:));

logicDemo=gcontains(allFiles,'DEMO_');

demoFiles=allFiles(logicDemo);
numDemoFiles=numel(demoFiles);

logicHelp=gcontains(allFiles,'HELP_');

helpFiles=allFiles(logicHelp);
numHelpFiles=numel(helpFiles);

logicStart=gcontains(allFiles,'GIBBON_product_page');

startFile=allFiles{logicStart};

%% ADJUST helptoc.XML
%TOP SECTION
T_top={};
T_top(1,1)={'<?xml version="1.0" encoding="utf-8"?>'};
T_top(2,1)={'<toc version="2.0">'};
T_top(3,1)={['    <tocitem target="',startFile,'">The GIBBON Toolbox']};

%GETTING STARTED SECTION
T_start=cell(3,1);
T_start{1}='        <tocitem target="GettingStarted.html" image="HelpIcon.GETTING_STARTED">Getting Started';
T_start{2}=['            <tocitem target="',startFile,'">',startFile(1:end-5),'</tocitem>',startFile(1:end-4)];
T_start{3}='        </tocitem>';

%FUNCTION HELP
T_help=cell(numHelpFiles+2,1);
T_help{1}='        <tocitem target="funclist.html" image="HelpIcon.FUNCTION">Functions';
for q=1:1:numHelpFiles
    currentFile=helpFiles{q};
    currentName=currentFile(1:end-5);
    T_help{q+1}=['            <tocitem target="',currentFile,'">',currentName(6:end),'</tocitem>',currentName];
end
T_help{end}='        </tocitem>';


%DEMO EXAMPLES
T_demo=cell(numDemoFiles+2,1);
T_demo{1}='        <tocitem target="gibbonExampes.html" image="HelpIcon.EXAMPLES">Examples';
for q=1:1:numDemoFiles
    currentFile=demoFiles{q};
    currentName=currentFile(1:end-5);
    T_demo{q+1}=['            <tocitem target="',currentFile,'">',currentName(6:end),'</tocitem>',currentName];
end
T_demo{end}='        </tocitem>';

%BOTTOM SECTION
T_bottom={};
T_bottom(1,1)={'    </tocitem>'};
T_bottom(2,1)={'</toc>'};

T=[T_top;T_start;T_help;T_demo;T_bottom];

saveName=fullfile(helpPath,'helptoc.xml');

cell2txtfile(saveName,T,0,0);

%% Add searchable help
addHelpSearch;

disp('Restart MATLAB to allow for the help and documentation integration changes to take effect');

 
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
