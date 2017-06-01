function createHelpDemoDocumentation

%% Get toolbox location

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
helpPath=fullfile(toolboxPath,'docs','html');

%% FIND HTML FILES

allFiles = dir(fullfile(helpPath,'*.html'));
allFiles={allFiles(1:end).name};
allFiles=sort(allFiles(:));

logicDemo=~cellfun(@isempty,strfind(allFiles,'DEMO_'));
demoFiles=allFiles(logicDemo);
numDemoFiles=numel(demoFiles);

logicHelp=~cellfun(@isempty,strfind(allFiles,'HELP_'));
helpFiles=allFiles(logicHelp);
numHelpFiles=numel(helpFiles);

indStart=find(~cellfun(@isempty,strfind(allFiles,'GIBBON_product_page')));
startFile=allFiles{indStart};

%% ADJUST helptoc.XML
%TOP SECTION
T_top={};
T_top(1,1)={'<?xml version="1.0" encoding="utf-8"?>'};
T_top(2,1)={'<toc version="2.0">'};
T_top(3,1)={'    <tocitem target="GIBBOM.html">The GIBBON Toolbox'};

%GETTING STARTED SECTION
T_start=cell(3,1);
T_start{1}=['        <tocitem target="GettingStarted.html" image="HelpIcon.GETTING_STARTED">Getting Started'];
T_start{2}=['            <tocitem target="',startFile,'">',startFile(1:end-5),'</tocitem>',startFile(1:end-4)];
T_start{3}=['        </tocitem>'];

%FUNCTION HELP
T_help=cell(numHelpFiles+2,1);
T_help{1}=['        <tocitem target="funclist.html" image="HelpIcon.FUNCTION">Functions'];
for q=1:1:numHelpFiles
    currentFile=helpFiles{q};
    currentName=currentFile(1:end-5);
    T_help{q+1}=['            <tocitem target="',currentFile,'">',currentName(6:end),'</tocitem>',currentName];
end
T_help{end}=['        </tocitem>'];


%DEMO EXAMPLES
T_demo=cell(numDemoFiles+2,1);
T_demo{1}=['        <tocitem target="gibbonExampes.html" image="HelpIcon.EXAMPLES">Examples'];
for q=1:1:numDemoFiles
    currentFile=demoFiles{q};
    currentName=currentFile(1:end-5);
    T_demo{q+1}=['            <tocitem target="',currentFile,'">',currentName(6:end),'</tocitem>',currentName];
end
T_demo{end}=['        </tocitem>'];

%BOTTOM SECTION
T_bottom={};
T_bottom(1,1)={'    </tocitem>'};
T_bottom(2,1)={'</toc>'};

T=[T_top;T_start;T_help;T_demo;T_bottom];

saveName=fullfile(helpPath,'helptoc.xml');

cell2txtfile(saveName,T,0);

%% Add searchable help
addHelpSearch;

disp('Restart MATLAB for the changes to take effect');

