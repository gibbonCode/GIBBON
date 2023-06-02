function installGibbon

% function installGibbon
%-------------------------------------------------------------------------
% Change log: 
% 2018/05/15 Added creation of temp folder if it does not exist
% 2020/11/24 Minor fix to suggested path names for FEBio and export_fig
%-------------------------------------------------------------------------

%% Add GIBBON library path so functions are known to use here
gibbonPath=fileparts(mfilename('fullpath')); %Get the GIBBON path
addpath(fullfile(gibbonPath,'lib')); %Add gibbon lib path so gibbon functions used here are known

%% Control parameters

W=60; %Scrollbar width
squareSizeMax=1200;

%% Set GIBBON settings

s=settings;

if hasGroup(s,'GIBBON') %Prior GIBBON settings group exists -> Check for setting
    if ~hasSetting(s.GIBBON,'vcw_profile') %Prior vcw_profile setting exists -> Update setting        
        addSetting(s.GIBBON,'vcw_profile'); % Add vcw_profile setting
        s.GIBBON.vcw_profile.PersonalValue='CAD'; %Make setting as personal setting
    end
else %No GIBBON settings at all -> Make group and setting
    addGroup(s,'GIBBON'); % Add GIBBON settings group
    addSetting(s.GIBBON,'vcw_profile'); % Add vcw_profile setting    
    s.GIBBON.vcw_profile.PersonalValue='CAD'; %Make setting as personal setting    
end

try  
    profileNameVCW=s.GIBBON.vcw_profile.PersonalValue; %Make setting as personal setting
catch
    profileNameVCW='CAD';
    s.GIBBON.vcw_profile.PersonalValue=profileNameVCW; %Make setting as personal setting
end

%% Open figure

close all; 

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


squareSize=round([min(screenSizeGroot(3:4)) min(screenSizeGroot(3:4))]); % width, height
squareSize(squareSize>squareSizeMax)=squareSizeMax;

hf.Units='pixels';
windowSize=[2*W (screenSizeGroot(4)-squareSize(2))-2*W squareSize(1) squareSize(2)]; % left bottom width height

hf.Position=windowSize; % left bottom width height

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

textPositionInstall=[W hf.Position(4)-W hf.Position(3)-W*2 round(W/1.5)];

top_title='Installing GIBBON';
hTextTitle = uicontrol(hf,'Style','text','String',top_title,...
    'Position',textPositionInstall,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',18,'FontWeight','normal');

hf.UserData.uihandles.hTextTitle=hTextTitle;

top_statement='Adding gibbon paths. Please wait...';
textPositionStatement=textPositionInstall;
textPositionStatement(2)=textPositionInstall(2)-W;

hTextStatement = uicontrol(hf,'Style','text','String',top_statement,...
    'Position',textPositionStatement,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','bold');

hf.UserData.uihandles.hTextStatement=hTextStatement;

%% Load images

imagePath=fullfile(gibbonPath,'docs','img');

M_overview=importdata(fullfile(imagePath,'GIBBON_overview.jpg'));
M_febio=importdata(fullfile(imagePath,'febio_banner.jpg'));
M_vcw=importdata(fullfile(imagePath,'viewControlConfiguration.jpg'));

hf.UserData.gibbonOverviewImage=M_overview;
hf.UserData.febioBanner=M_febio;
hf.UserData.viewControlImage=M_vcw;

%% Prepare axis and add GIBBON overview image

figure(hf);
hImage=image(hf.UserData.gibbonOverviewImage);
axis tight; axis equal; axis off;
hAxis=gca;
hAxis.Units='pixels';

% [left bottom width height]
axisHeight=hf.Position(4)-7*W;
axisWidth=hf.Position(3)-2*W;
hAxis.Position=[round((hf.Position(3)-axisWidth)/2) round((hf.Position(4)-axisHeight))-6*W axisWidth axisHeight];

hf.UserData.uihandles.hAxis=hAxis;
hf.UserData.hImage=hImage;

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

for q=1:1:numel(pathNames)
     pathNameNow=pathNames{q};     
     addpath(pathNameNow); %Add path     
     hTextStatement.String=[top_statement,' ',num2str(round(100*q/numel(pathNames))),'% complete'];
     drawnow;
end

hTextStatement.String='Done adding toolbox paths';
drawnow;
pause(0.5);

%% Add FEBio path

delete(hf.UserData.hImage);
hf.UserData.hImage=image(hf.UserData.febioBanner);
axis tight; axis equal; axis off;

hTextStatement.String='Please provide the full path to the FEBio executable. Leave blank if FEBio is not needed.';

if ispc
    hTextInfoStringDefault='On Windows likely similar to: C:\Program Files\FEBioStudio2\febio\FEBio4.exe';
elseif ismac
    hTextInfoStringDefault='On MacOS likely similar to: /Applications/FEBioStudio/FEBioStudio.app/Contents/MacOS/febio4';
else
    hTextInfoStringDefault='On Linux likely similar to: /home/<UserName>/FEBioStudio/bin/febio4';
end

hTextInfo1 = uicontrol(hf,'Style','text','String',hTextInfoStringDefault,...
    'Position',[W hf.Position(4)-W*3 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','normal');

FEBioPath=getFEBioPath;

hTextInput1 = uicontrol(hf,'Style','edit','String',FEBioPath,...
    'Position',[W hf.Position(4)-W*3.5 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12);

hf.UserData.uihandles.hTextInfo1=hTextInfo1;
hf.UserData.uihandles.hTextInput1=hTextInput1;

%% Create push button
hf.UserData.pathDefinitionsDone=0;

hBrowseButton = uicontrol('Style', 'pushbutton', 'String', 'Browse','Position',[W hf.Position(4)-W*4.5   5*W round(W/1.5)],'Callback',{@getFEBioExecPath,{hf}},'FontSize',12);
hf.UserData.uihandles.hBrowseButton=hBrowseButton;

hconfirmButton = uicontrol('Style', 'pushbutton', 'String', 'Confirm path','Position',[W hf.Position(4)-W*5.5   5*W round(W/1.5)],'Callback',{@setThirdpartyPaths,{hf}},'FontSize',12);
hf.UserData.uihandles.hconfirmButton=hconfirmButton;

% Wait for path definitions to be set

while hf.UserData.pathDefinitionsDone==0
    pause(0.1);
end

%% Add image

hf.UserData.VCWOptionDone=0;

indCurrent=find(strcmp(profileNameVCW,{'CAD','febio','touchpad'}));

hf.UserData.uihandles.hTextStatement.String=['Set+test view control widget button mapping profile. Current setting: ',profileNameVCW];

delete(hf.UserData.hImage);
hf.UserData.hImage=image(hf.UserData.viewControlImage);
axis tight; axis equal; axis off;

hVCWSelect = uicontrol('Style', 'popupmenu', 'String', {'CAD','febio','touchpad'},'Value',indCurrent,'Position',[W hf.Position(4)-W*3.5   5*W round(W/1.5)],'Callback',{@setVCWOption,{hf}},'FontSize',12);
hf.UserData.uihandles.hVCWSelect=hVCWSelect;

hconfirmButton = uicontrol('Style', 'pushbutton', 'String', 'Confirm view profile','Position',[W hf.Position(4)-W*4.5   5*W round(W/1.5)],'Callback',{@VCWOptionDone,{hf}},'FontSize',12);
hf.UserData.uihandles.hconfirmButton=hconfirmButton;

createTestFigure(hf,profileNameVCW);

while hf.UserData.VCWOptionDone==0
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

%%

hf.UserData.uihandles.hTextTitle.String='Finished! GIBBON is installed.';

hf.UserData.uihandles.hTextStatement.String='Use gdoc to open the GIBBON help from the command window, or visit https://www.gibboncode.org/Documentation/';

hDoneButton = uicontrol('Style', 'pushbutton', 'String', 'Done!','Position',[W hf.Position(4)-W*3.5   5*W round(W/1.5)],'Callback',{@closeDone,{hf}},'FontSize',12);
hf.UserData.uihandles.hBrowseButton=hDoneButton;

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

if isfield(hf.UserData.uihandles,'hBrowseButton')
    hf.UserData.uihandles.hBrowseButton.Position=     [W hf.Position(4)-W*4.5   5*W                round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hconfirmButton')
    hf.UserData.uihandles.hconfirmButton.Position=     [W hf.Position(4)-W*5.5   5*W                round(W/1.5)];
end

if isfield(hf.UserData.uihandles,'hDoneButton')
    hf.UserData.uihandles.hDoneButton.Position=     [W hf.Position(4)-W*3.5   5*W                round(W/1.5)];
end

end

%%

function setThirdpartyPaths(~,~,inputCell)
hf=inputCell{1};

t=hf.UserData.uihandles.hTextInput1.String;
if ~isempty(t)
    setFEBioPath(t); %Set FEBio path in config file
end

%Clean-up guide components
delete(hf.UserData.uihandles.hTextInfo1)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInfo1');
delete(hf.UserData.uihandles.hTextInput1)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hTextInput1');

delete(hf.UserData.uihandles.hBrowseButton)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hBrowseButton');

delete(hf.UserData.uihandles.hconfirmButton)
hf.UserData.uihandles=rmfield(hf.UserData.uihandles,'hconfirmButton');

hf.UserData.uihandles.hTextStatement.String='Done, adding third party paths';
drawnow;

% Saving path definitions
hf.UserData.uihandles.hTextStatement.String='Saving path definitions';
drawnow;

savepath;

hf.UserData.uihandles.hTextStatement.String='Done, saving path definitions';
drawnow;

hf.UserData.pathDefinitionsDone=1;

end

%% Browse to exec file

function getFEBioExecPath(~,~,inputCell)

hf=inputCell{1};
[fileName,filePath]=uigetfile('*','Select the FEBio executable file');
if fileName==0 %Selection was cancelled
    hf.UserData.uihandles.hTextInput1.String='';
else %Selection was made
    hf.UserData.uihandles.hTextInput1.String=fullfile(filePath,fileName);
end
drawnow; 

end


%% Set VCW option

function setVCWOption(~,~,inputCell)

hf=inputCell{1};

s=hf.UserData.uihandles.hVCWSelect.String;
profileNameVCW=s{hf.UserData.uihandles.hVCWSelect.Value};
setViewProfile(profileNameVCW);

hf.UserData.uihandles.hTextStatement.String=['Set+test view control widget button mapping profile. Current setting: ',profileNameVCW];

createTestFigure(hf,profileNameVCW); 

end

function VCWOptionDone(~,~,inputCell)

hf=inputCell{1};

hf.UserData.uihandles.hTextStatement.String='Done, setting view control profile';
drawnow;

hf.UserData.VCWOptionDone=1;

%Clean-up guide components
hf.UserData.uihandles.hTextStatement.String='';
delete(hf.UserData.uihandles.hVCWSelect); %Delete toggle input
delete(hf.UserData.uihandles.hconfirmButton); %Delete confirm button

if isfield(hf.UserData,'hf2')
    close(hf.UserData.hf2);
end

delete(hf.UserData.hImage);
hf.UserData.hImage=image(hf.UserData.gibbonOverviewImage);
axis tight; axis equal; axis off;

end

function createTestFigure(hf,profileNameVCW)

if isfield(hf.UserData,'hf2')
    close(hf.UserData.hf2);
end

gibbonPath=fileparts(mfilename('fullpath')); %Get the GIBBON path
fileName=fullfile(gibbonPath,'data','OBJ','gibbon.obj');

objStruct=import_obj(fileName);
F=objStruct.F;
V=objStruct.V;
C=objStruct.C;

hf.UserData.hf2=figure; hold on; 
s=['Testing VCW profile: ',profileNameVCW];
title(s,'FontSize',25); hf.UserData.hf2.Name=s;
vcw(hf.UserData.hf2,profileNameVCW);
hp=gpatch(F,V,C,'none'); 
axisGeom; camlight headlight; gdrawnow; 

end

%% Close/done

function closeDone(~,~,inputCell)

hf=inputCell{1};
close(hf);

try
    gdoc
catch 
    
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

