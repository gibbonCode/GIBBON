function installGibbon

%% Add GIBBON library path so functions are known to use here
gibbonPath=fileparts(mfilename('fullpath')); %Get the GIBBON path
addpath(fullfile(gibbonPath,'lib')); %Add gibbon lib path so gibbon functions used here are known

%% Settings

W=40; %Scrollbar width

%% Open figure

figStruct.Name='Installing GIBBON';


hf = cFigure(figStruct);
hf.NumberTitle='off';

%Remove tool and menu bars
ht = findobj(allchild(hf),'flat','Type','uitoolbar');
delete(ht);
ht = findobj(allchild(hf),'flat','Type','uimenu');
delete(ht);

set(hf,'ResizeFcn',{@resizeAll,{hf}});

% [left bottom width height]

hf.UserData.W=W;

%% Initilize text fields

top_title='Installing GIBBON';
hTextTitle = uicontrol(hf,'Style','text','String',top_title,...
    'Position',[W hf.Position(4)-W hf.Position(3)-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',20,'FontWeight','normal');

hf.UserData.uihandles.hTextTitle=hTextTitle;

top_statement='Adding gibbon paths. Please wait...';
hTextStatement = uicontrol(hf,'Style','text','String',top_statement,...
    'Position',[W hf.Position(4)-2*W hf.Position(3)-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',16,'FontWeight','bold');

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

%% Adding paths

[pathNames]=getSubPaths(gibbonPath);

%Remove "hidden" folders with . in path name
logicKeep=~contains(pathNames,'.');
pathNames=pathNames(logicKeep);

hw = waitbar(0,'Adding toolbox paths');
for q=1:1:numel(pathNames)
     pathNameNow=pathNames{q};     
     addpath(pathNameNow); %Add path     
     waitbar(q/numel(pathNames),hw);
end
close(hw)

hTextStatement.String='Done adding toolbox paths';

%% Add 3rd party paths

% % disp('-> Adding 3rd party paths');
% %
% % prompt = {'FEBio full path to program (e.g. .../bin/FEBio2.lnx64 or ...\bin\FEBio2.exe):','export_fig path:'};
% % dlg_title = 'Path definitions (leave empty if not used)';
% %
% % FEBioPath=getFEBioPath;
% % % if isempty(FEBioPath)
% % %     if ispc
% % %         FEBioPath='';
% % %     elseif isunix
% % %         FEBioPath='';
% % %     end
% % % end
% %
% % exportFigPath=fileparts(which('export_fig'));
% % defaultOptions = {FEBioPath,exportFigPath};
% %
% % s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
% %
% % Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
% %
% % if ~isempty(Q)
% %     if ~isempty(Q{1})
% %         setFEBioPath(Q{1}); %Set FEBio path in config file
% %     end
% %     if ~isempty(Q{2})
% %         addpath(Q{2}); %Add export_fig to the path
% %     end
% % end
%
%
% %%

hTextStatement.String='Please provide 3rd party package locations.';

hTextInfoStringDefault='Full path to FEBio excutable (leave blank if not needed):';

hTextInfo1 = uicontrol(hf,'Style','text','String',hTextInfoStringDefault,...
    'Position',[W hf.Position(4)-W*3 round(hf.Position(3))-W*2 round(W/1.5)],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','bold');

FEBioPath=getFEBioPath;

if isempty(FEBioPath)
    if ispc
        FEBioPath='e.g. C:\Program Files\febio-2.6.4\bin\FEBio2.exe';
    else
        FEBioPath='e.g. /home/febio-2.6.4/bin/febio2.lnx64';
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
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left','FontSize',12,'FontWeight','bold');

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

%% Create push button
hconfirmButton = uicontrol('Style', 'pushbutton', 'String', 'Set paths','Position',[W hf.Position(4)-W*6 5*W round(W/1.5)],'Callback',{@setThirdpartyPaths,{hf}},'FontSize',12);
hf.UserData.uihandles.hconfirmButton=hconfirmButton;

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
pause(0.1);

%% Saving path definitions
hf.UserData.uihandles.hTextStatement.String='Saving path definitions';
savepath;
hf.UserData.uihandles.hTextStatement.String='Done, saving path definitions';
pause(0.1);

%% Integrating help/documentations

hf.UserData.uihandles.hTextStatement.String='Integrating help';
createHelpDemoDocumentation;

hf.UserData.uihandles.hTextStatement.String='Restart MATLAB to allow for the help and documentation integration changes to take effect';

hf.UserData.uihandles.hTextTitle.String='Finished GIBBON is installed. Feel free to close this window';

end
