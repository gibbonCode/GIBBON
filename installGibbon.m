function installGibbon

% function installGibbon
%-------------------------------------------------------------------------
% Change log: 
% 2018/05/15 Added creation of temp folder if it does not exist
%-------------------------------------------------------------------------

%% Add GIBBON library path so functions are known to use here
gibbonPath=fileparts(mfilename('fullpath')); %Get the GIBBON path
addpath(fullfile(gibbonPath,'lib')); %Add gibbon lib path so gibbon functions used here are known

%% Settings

W=60; %Scrollbar width

%% Open figure

%Force groot units to be pixels
graphicalRoot=groot;
grootUnits=graphicalRoot.Units;
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units='pixels';
end

screenSizeGroot = get(groot,'ScreenSize');

figStruct.Name='Installing GIBBON';

hf = cFigure(figStruct);
hf.NumberTitle='off';

figSize=round([min(screenSizeGroot(3:4)) min(screenSizeGroot(3:4))*0.9]); % width, height

hf.Units='pixels';
hf.Position=[(screenSizeGroot(3)-figSize(1))/2 (screenSizeGroot(4)-figSize(2))/2  figSize(1) figSize(2)]; % left bottom width height

%Remove tool and menu bars
ht = findobj(allchild(hf),'flat','Type','uitoolbar');
delete(ht);
ht = findobj(allchild(hf),'flat','Type','uimenu');
delete(ht);

set(hf,'ResizeFcn',{@resizeAll,{hf}});

hf.UserData.W=W;
hf.UserData.gibbonPath=gibbonPath;

% Reset groot units if a change was needed
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units=grootUnits;
end

%% Initilize text fields

top_title='Installing GIBBON';
hTextTitle = uicontrol(hf,'Style','text','String',top_title,...
    'Position',[W hf.Position(4)-W hf.Position(3)-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',18,'FontWeight','normal');

hf.UserData.uihandles.hTextTitle=hTextTitle;

top_statement='Adding gibbon paths. Please wait...';
hTextStatement = uicontrol(hf,'Style','text','String',top_statement,...
    'Position',[W hf.Position(4)-2*W hf.Position(3)-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','bold');

hf.UserData.uihandles.hTextStatement=hTextStatement;

%% Add image

imagePath=fullfile(gibbonPath,'docs','img','GIBBON_overview.jpg');
M=importdata(imagePath);

figure(hf);
image(M);
axis tight; axis equal; axis off;
hAxis=gca;

hAxis.Units='pixels';
% [left bottom width height]
axisHeight=hf.Position(4)-7*W;
axisWidth=hf.Position(3)-2*W;
hAxis.Position=[round((hf.Position(3)-axisWidth)/2) round((hf.Position(4)-axisHeight))-6*W axisWidth axisHeight];

hf.UserData.uihandles.hAxis=hAxis;

%%
drawnow;

%% Adding temp folder if it does not exist

tempPath=fullfile(gibbonPath,'data','temp');

if ~exist(tempPath,'file')
    mkdir(tempPath)
end

%% Adding paths

[pathNames]=getSubPaths(gibbonPath);

%Remove "hidden" folders with . in path name
logicKeep=~gcontains(pathNames,'.');
pathNames=pathNames(logicKeep);

hw = waitbar(0,'Adding toolbox paths');
for q=1:1:numel(pathNames)
     pathNameNow=pathNames{q};     
     addpath(pathNameNow); %Add path     
     waitbar(q/numel(pathNames),hw);
end
close(hw)

hTextStatement.String='Done adding toolbox paths';
drawnow;

%% Add 3rd party paths
hTextStatement.String='Please provide 3rd party package locations. Leave blank if associated features are not needed.';

hTextInfoStringDefault='Full path to FEBio excutable, e.g. /home/febio-2.7.1/bin/febio2.lnx64 (or FEBio2.exe for windows):';

hTextInfo1 = uicontrol(hf,'Style','text','String',hTextInfoStringDefault,...
    'Position',[W hf.Position(4)-W*3 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','normal');

FEBioPath=getFEBioPath;

if isempty(FEBioPath)
    if ispc
        FEBioPath='C:\Program Files\febio-2.8.3\bin\FEBio2.exe';
    else
        FEBioPath='/home/febio-2.8.3/bin/febio2.lnx64';
    end
end

hTextInput1 = uicontrol(hf,'Style','edit','String',FEBioPath,...
    'Position',[W hf.Position(4)-W*3.5 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12);

hf.UserData.uihandles.hTextInfo1=hTextInfo1;
hf.UserData.uihandles.hTextInput1=hTextInput1;

%%
hTextInfoStringDefault='Full path to export_fig:';

hTextInfo2 = uicontrol(hf,'Style','text','String',hTextInfoStringDefault,...
    'Position',[W hf.Position(4)-W*4.5 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','normal');

exportFigPath=fileparts(which('export_fig'));

if isempty(FEBioPath)
    if ispc
        exportFigPath='e.g. ...\MATLAB\export_fig';
    else
        exportFigPath='e.g. ...C:\Program Files\febio-2.5.2\bin';
    end
end

hTextInput2 = uicontrol(hf,'Style','edit','String',exportFigPath,...
    'Position',[W hf.Position(4)-W*5 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12);

hf.UserData.uihandles.hTextInfo2=hTextInfo2;
hf.UserData.uihandles.hTextInput2=hTextInput2;

%% 

hTextStatement.String='Please provide 3rd party package locations.';

%% Create push button
hf.UserData.pathDefinitionsDone=0;

hconfirmButton = uicontrol('Style', 'pushbutton', 'String', 'Set paths','Position',[W hf.Position(4)-W*6 5*W round(W/1.5)],'Callback',{@setThirdpartyPaths,{hf}},'FontSize',12);
hf.UserData.uihandles.hconfirmButton=hconfirmButton;

%% Wait for path definitions to be set

while hf.UserData.pathDefinitionsDone==0
    pause(0.1);
end

%% Unzip compressed data
hf.UserData.uihandles.hTextStatement.String='Unzipping compressed data'; 
drawnow;

dataFolder=fullfile(hf.UserData.gibbonPath,'data');
unzipAll(dataFolder,1);

hf.UserData.uihandles.hTextStatement.String='Done unzipping data'; drawnow;
pause(0.5); 

%% Integrating help/documentations
hf.UserData.uihandles.hTextStatement.String='Integrating help'; drawnow;
createHelpDemoDocumentation;
hf.UserData.uihandles.hTextStatement.String='Done integrating help'; drawnow;
pause(0.5); 

%% Compiling MEX files
% hf.UserData.uihandles.hTextStatement.String='Compiling mex files (required for faster geodesic remeshing)'; drawnow;
% 
% answer = questdlg('Would you like to try compile mex files now?','Compiling mex files','Yes','No','Yes');
% 
% switch answer
%     case 'Yes'
%         compile_mex;
%         hf.UserData.uihandles.hTextStatement.String='Done compiling mex files'; drawnow;
%     otherwise
%         hf.UserData.uihandles.hTextStatement.String='Skipped compiling mex files'; drawnow;     
% end
% pause(0.5); 

%%

hf.UserData.uihandles.hTextStatement.String='Restart MATLAB to allow for the help and documentation integration changes to take effect';
hf.UserData.uihandles.hTextTitle.String='Finished! GIBBON is installed. Feel free to close this window';

end

%% Resizing

function resizeAll(~,~,inputCell)
hf=inputCell{1};
W=hf.UserData.W;

%Reposition axis
axisHeight=hf.Position(4)-7*W;
axisWidth=hf.Position(3)-2*W;
hf.UserData.uihandles.hAxis.Position=[round((hf.Position(3)-axisWidth)/2) round((hf.Position(4)-axisHeight))-6*W axisWidth axisHeight];

%Reposition gui items
if isfield(hf.UserData.uihandles,'hTextTitle')
    hf.UserData.uihandles.hTextTitle.Position=      [W hf.Position(4)-W     hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hTextStatement')
    hf.UserData.uihandles.hTextStatement.Position=  [W hf.Position(4)-2*W   hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hTextInfo1')
    hf.UserData.uihandles.hTextInfo1.Position=      [W hf.Position(4)-3*W hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hTextInput1')
    hf.UserData.uihandles.hTextInput1.Position=     [W hf.Position(4)-3.5*W   hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hTextInfo2')
    hf.UserData.uihandles.hTextInfo2.Position=      [W hf.Position(4)-4.5*W hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hTextInput2')
    hf.UserData.uihandles.hTextInput2.Position=     [W hf.Position(4)-5*W   hf.Position(3)-W*2 round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hconfirmButton')
    hf.UserData.uihandles.hconfirmButton.Position=     [W hf.Position(4)-W*6   5*W                round(W/1.5)];
end


end

%%

function setThirdpartyPaths(~,~,inputCell)
hf=inputCell{1};

t=hf.UserData.uihandles.hTextInput1.String;
if ~isempty(t)
    setFEBioPath(t); %Set FEBio path in config file
end

t=hf.UserData.uihandles.hTextInput2.String;
if ~isempty(t)
    addpath(t); %Add export_fig to the path
end

delete(hf.UserData.uihandles.hTextInfo1)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInfo1');
delete(hf.UserData.uihandles.hTextInput1)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInput1');

delete(hf.UserData.uihandles.hTextInfo2)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInfo2');
delete(hf.UserData.uihandles.hTextInput2)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInput2');

delete(hf.UserData.uihandles.hconfirmButton)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hconfirmButton');

hf.UserData.uihandles.hTextStatement.String='Done, adding third party paths';
drawnow;

%% Saving path definitions
hf.UserData.uihandles.hTextStatement.String='Saving path definitions';
drawnow;

savepath;

hf.UserData.uihandles.hTextStatement.String='Done, saving path definitions';
drawnow;

hf.UserData.pathDefinitionsDone=1;

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

