function uninstallGibbon

% function uninstallGibbon
%-------------------------------------------------------------------------
% Change log: 
% 2022/03/16 Created based on installGibbon
%-------------------------------------------------------------------------

%% Add GIBBON library path so functions are known to use here
gibbonPath=fileparts(mfilename('fullpath')); %Get the GIBBON path
addpath(fullfile(gibbonPath,'lib')); %Add gibbon lib path so gibbon functions used here are known

%% Settings

W=60; %Scrollbar width
squareSizeMax=1200;

%% Open figure

%Force groot units to be pixels
graphicalRoot=groot;
grootUnits=graphicalRoot.Units;
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units='pixels';
end

screenSizeGroot = get(groot,'ScreenSize');

figStruct.Name='uninstalling GIBBON';

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

top_title='uninstalling GIBBON';
hTextTitle = uicontrol(hf,'Style','text','String',top_title,...
    'Position',textPositionInstall,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',18,'FontWeight','normal');

hf.UserData.uihandles.hTextTitle=hTextTitle;

top_statement='Removing gibbon paths. Please wait...';
textPositionStatement=textPositionInstall;
textPositionStatement(2)=textPositionInstall(2)-W;

hTextStatement = uicontrol(hf,'Style','text','String',top_statement,...
    'Position',textPositionStatement,...
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

if exist(tempPath,'file')
    rmdir(tempPath,'s');
end

%% Adding paths

[pathNames]=getSubPaths(gibbonPath);

%Remove "hidden" folders with . in path name
logicKeep=~gcontains(pathNames,'.');
pathNames=pathNames(logicKeep);

for q=1:1:numel(pathNames)
     pathNameNow=pathNames{q};     
     rmpath(pathNameNow); %Remove path
     hTextStatement.String=[top_statement,' ',num2str(round(100*q/numel(pathNames))),'% complete'];
     drawnow;
end
rmpath(fullfile(gibbonPath,'lib'));

hTextStatement.String='Done removing toolbox paths';
drawnow;
pause(0.5);

%% Unzip compressed data
hf.UserData.uihandles.hTextStatement.String='Removing unzipped data'; 
drawnow;

zipFolder=fullfile(hf.UserData.gibbonPath,'data','DICOM/0001_human_calf');

files = dir(fullfile(zipFolder,'*.dcm'));
files={files(1:end).name};
files=sort(files(:));
if ~isempty(files)
    for q=1:numel(files)
        delete(fullfile(zipFolder,files{q}));
    end
end

hf.UserData.uihandles.hTextStatement.String='Done removing unzipped data'; drawnow;
pause(0.5); 

%% Removing documentation integrating

% hf.UserData.uihandles.hTextStatement.String='Removing documentation'; drawnow;
% 
% hf.UserData.uihandles.hTextStatement.String='Done removing documentation'; drawnow;
% pause(0.5); 

%%

hf.UserData.uihandles.hTextTitle.String='Done. GIBBON has been uninstalled.';

hf.UserData.uihandles.hTextStatement.String='Documentation integration removal is not supported yet. ';

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


%% Close/done

function closeDone(~,~,inputCell)

hf=inputCell{1};
close(hf);

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

