function hf=imx(varargin)

%%


%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
    case 2
        M=varargin{1};
        v=varargin{2};
end

M=double(M);

siz=size(M);

%Cope with 2D image sets
if ismatrix(M)
    siz(3)=1;
end

%%

%Get image coordinates
[J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));

%Convert to cartesian coordinates using voxel size if provided
[X,Y,Z]=im2cart(I,J,K,v);

%% Plot settings
fontColor='w';
fontSize=20;
cMap=gray(250);

figStruct.Name='Image Segmentation Widget'; %Figure name
figStruct.Color=0.*ones(1,3); %Figure background color
figStruct.ColorDef='black'; %Setting colordefinitions to black
% figStruct.ScreenOffset=100; %Setting spacing of figure with respect to screen edges

%% Defining row, column and slice indicices for slice patching
sliceIndexI=round(siz(1)/2); %(close to) middle row
sliceIndexJ=round(siz(2)/2); %(close to) middle column
sliceIndexK=round(siz(3)/2); %(close to) middle slice

Tf=[0 100]; %Threshold

nTickMajor=20;
tickSizeMajor_I=round(siz(1)/nTickMajor);
tickSizeMajor_J=round(siz(2)/nTickMajor);
tickSizeMajor_K=round(siz(3)/nTickMajor);

%% Initialize display

[ax,ay,az]=im2cart([siz(1)+1 0],[siz(2)+1 0],[siz(3)+1 0],v);

hf=cFigure(figStruct);
axis equal; axis tight; axis([ax(2) ax(1) ay(2) ay(1) az(2) az(1)]); view(3);  axis vis3d; grid on; box on; hold on;
view(0,90);
hAxis=gca;
axis_units=hAxis.Units;
hAxis.Units='pixels';
hAxis.Position(2)=50;
hAxis.Position(4)=hf.Position(4)-100;
hAxis.Units=axis_units;

colormap(cMap);
hColorBar=colorbar;
hColorBar.Units='pixels';
hColorBar.AxisLocation='in';
hColorBar.Location='south';

hColorBar.Position(1)=4*60;
hColorBar.Position(2)=10;
hColorBar.Position(3)=hf.Position(3)-(4*60)-10;
hColorBar.Position(4)=25;
hColorBar.Label.String = '';
% hColorBar.FontSize=10;
set(gca,'fontSize',fontSize);
%
% hColorBar=colorbar;
% hColorBar.Units='pixels';
% hColorBar.AxisLocation='in';
% % hColorBar.Position(1)=hf.Position(3)-50;
% % hColorBar.Position(2)=50;
% % hColorBar.Position(3)=50/2;
% % hColorBar.Position(4)=hf.Position(4)-100;

caxis([min(M(:)) max(M(:))]);

set(gca,'fontSize',fontSize);
drawnow;

w=50; %Scrollbar width
jSlider_I = javax.swing.JSlider(1,siz(1));
javacomponent(jSlider_I,[0,0,w,round(hf.Position(4))]);
set(jSlider_I, 'MajorTickSpacing',tickSizeMajor_I, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_I,1}},'Orientation',jSlider_I.VERTICAL);

jSlider_J = javax.swing.JSlider(1,siz(2));
javacomponent(jSlider_J,[1*w,0,w,round(hf.Position(4))]);
set(jSlider_J, 'MajorTickSpacing',tickSizeMajor_J, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_J,2}},'Orientation',jSlider_J.VERTICAL);

jSlider_K = javax.swing.JSlider(1,siz(3));
javacomponent(jSlider_K,[2*w,0,w,round(hf.Position(4))]);
set(jSlider_K, 'MajorTickSpacing',tickSizeMajor_K, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_K,3}},'Orientation',jSlider_K.VERTICAL);

jSlider_T = com.jidesoft.swing.RangeSlider(0,100,Tf(1),Tf(2));  % min,max,low,high
javacomponent(jSlider_T,[3*w,0,w,round(hf.Position(4))]);
set(jSlider_T, 'MajorTickSpacing',25, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',false, 'StateChangedCallback',{@setThresholdFunc,{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K}},'Orientation',jSlider_T.VERTICAL);

%% Set resize function

set(hf,'ResizeFcn',{@setScrollSizeFunc,{hf,w,jSlider_T,jSlider_I,jSlider_J,jSlider_K}});

%% Initialize figure callbacks
% set(hf, 'WindowButtonDownFcn', {@FigMouseDown,hf}, ...
%     'WindowButtonUpFcn', {@FigMouseUp,hf}, ...
%     'KeyPressFcn', {@FigKeyPress,hf}, ...
%     'WindowScrollWheelFcn', {@FigScroll,hf}, ...
%     'BusyAction', 'cancel');

set(hf,'KeyPressFcn', {@figKeyPressFunc,{hf}},'WindowButtonDownFcn', {@figMouseDown,{hf}},'WindowButtonUpFcn', {@mouseup,hf});

%% Initialise buttons

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
iconPath=fullfile(toolboxPath,'icons');

hb = findall(hf,'Type','uitoolbar');
% hb = uitoolbar(hf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help button

%get icon
D=importdata(fullfile(iconPath,'help.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Help','CData',S,'Tag','help_button','ClickedCallback',@helpFunc,'Separator','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save button

%get icon
D=importdata(fullfile(iconPath,'save.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Save','CData',S,'Tag','save_button','ClickedCallback',{@saveFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load button

%get icon
D=importdata(fullfile(iconPath,'load.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Load','CData',S,'Tag','load_button','ClickedCallback',{@loadFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Colorbar button

%get icon
D=importdata(fullfile(iconPath,'colorbar.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','colorbar','CData',S,'Tag','colorbar_button','ClickedCallback',@colorbarFunc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ghost button

%get icon
D=importdata(fullfile(iconPath,'ghost.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Transparency','CData',S,'Tag','ghost_button','ClickedCallback',{@ghostFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%showHide button

%get icon
D=importdata(fullfile(iconPath,'eye.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','showHide','CData',S,'Tag','showHide_button','ClickedCallback',{@showHideFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Home button

%get icon
D=importdata(fullfile(iconPath,'home.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Home','CData',S,'Tag','home_button','ClickedCallback',{@homeFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample button

%get icon
D=importdata(fullfile(iconPath,'sample.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hSample=uitoggletool(hb,'TooltipString','Sample contour','CData',S,'Tag','sample_button','ClickedCallback',{@sampleFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cut button

%get icon
D=importdata(fullfile(iconPath,'cut.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hCut=uitoggletool(hb,'TooltipString','Cut contours','CData',S,'Tag','cut_button','ClickedCallback',{@cutFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw button

%get icon
D=importdata(fullfile(iconPath,'draw.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hDraw=uitoggletool(hb,'TooltipString','Draw manual contour','CData',S,'Tag','draw_button','ClickedCallback',{@drawFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select button

%get icon
D=importdata(fullfile(iconPath,'select.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hSelect=uitoggletool(hb,'TooltipString','Select','CData',S,'Tag','select_button','ClickedCallback',{@selectFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Delete contour button

%get icon
D=importdata(fullfile(iconPath,'delete.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hDelete=uitoggletool(hb,'TooltipString','Delete','CData',S,'Tag','delete_button','ClickedCallback',{@deleteFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reset contours button

%get icon
D=importdata(fullfile(iconPath,'reset.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Reset','CData',S,'Tag','reset_button','ClickedCallback',{@resetFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert button

%get icon
D=importdata(fullfile(iconPath,'convert.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
% S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Convert','CData',S,'Tag','convert_button','ClickedCallback',{@convertFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Close contour button

%get icon
D=importdata(fullfile(iconPath,'closePolygon.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Close','CData',S,'Tag','close_button','ClickedCallback',{@closeFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smooth button

%get icon
D=importdata(fullfile(iconPath,'smooth.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hSmooth=uitoggletool(hb,'TooltipString','Smooth','CData',S,'Tag','smooth_button','ClickedCallback',{@smoothFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ellipse button

%get icon
D=importdata(fullfile(iconPath,'ellipse.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
% S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hEllipse=uitoggletool(hb,'TooltipString','Ellipse','CData',S,'Tag','ellipse_button','ClickedCallback',{@ellipseFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grow button

%get icon
D=importdata(fullfile(iconPath,'grow.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hGrow=uitoggletool(hb,'TooltipString','Grow','CData',S,'Tag','grow_button','ClickedCallback',{@growShrinkFunc,{hf,1}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shrink button

%get icon
D=importdata(fullfile(iconPath,'shrink.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
hShrink=uitoggletool(hb,'TooltipString','Shrink','CData',S,'Tag','shrink_button','ClickedCallback',{@growShrinkFunc,{hf,-1}});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Target button
%
% %get icon
% D=importdata(fullfile(iconPath,'target.jpg'));
% S=double(D);
% S=S-min(S(:));
% S=S./max(S(:));
% S(S==1)=NaN;
% if size(S,3)==1
%     S=repmat(S,[1 1 3]);
% end
% % Create a uipushtool in the toolbar
% uitoggletool(hb,'TooltipString','Target','CData',S,'Tag','target_button',);%,'ClickedCallback',{@loadFunc,{hf}});

%% Popup
% hPopUp1 = uicontrol('Style', 'popup',...
%     'String', {'gray','gjet','jet','hsv'},...
%     'Position', [hf.Position(3)-100 hf.Position(4)-20 100 20],...
%     'Callback', @funcColormap);

% %% Text overlay
%
% topTextString='Vertical Exaggeration';
% pos=hPopUp1.Position; pos(2)=pos(2)-pos(4);
% hTopText = uicontrol('Style','text','Position',pos,'String',topTextString,'BackgroundColor',0.5*ones(1,3),'FontSize',10,'ForegroundColor',abs(1-figStruct.Color));

%% Set figure UserData

hf.UserData.M=M;
hf.UserData.Name=figStruct.Name;
hf.UserData.v=v;

hf.UserData.X=X;
hf.UserData.Y=Y;
hf.UserData.Z=Z;

hf.UserData.patchTypes={'si','sj','sk'};
hf.UserData.sliceIndices=[sliceIndexI sliceIndexJ sliceIndexK];
hf.UserData.contourSlice=sliceIndexK;
hf.UserData.fontColor=fontColor;
hf.UserData.faceAlpha=0.8;
hf.UserData.lineColors=gjet(4);

hf.UserData.ButtonHandles.Sample=hSample;
hf.UserData.ButtonHandles.Cut=hCut;
hf.UserData.ButtonHandles.Draw=hDraw;
hf.UserData.ButtonHandles.Select=hSelect;
hf.UserData.ButtonHandles.Delete=hDelete;
hf.UserData.ButtonHandles.hSmooth=hSmooth;
hf.UserData.ButtonHandles.hEllipse=hEllipse;
hf.UserData.ButtonHandles.hGrow=hGrow;
hf.UserData.ButtonHandles.hShrink=hShrink;

hf.UserData.colorBarhandle=hColorBar;
% hf.UserData.hPopUp1=hPopUp1;
hf.UserData.hAxis=hAxis;

hf.UserData.sketchContourHandle=[];
hf.UserData.contourSetHandle=[];

hf.UserData.sliderHandles={jSlider_I,jSlider_J,jSlider_K,jSlider_T};

hf.UserData.savePath=fullfile(cd,'imseg'); %Save path
hf.UserData.saveName='imseg'; %Save name

hf.UserData.sketchContour={};
hf.UserData.ContourSet=repmat({{[]}},1,siz(3));

cMapContours=gjet(4);
hf.UserData.colormapSketch=[cMapContours(3,:); cMapContours(4,:)]; %red to yellow
hf.UserData.colormapSet=[cMapContours(1,:); cMapContours(2,:)]; %blue to green

% Store current settings
hf.UserData.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.BusyAction=hf.BusyAction;

hf.UserData.csapsSmoothPar=0.5; 

hf.UserData.showAll=1; 

%% Initialize slider locations
set(jSlider_I,'Value',sliceIndexI);
set(jSlider_J,'Value',sliceIndexJ);
set(jSlider_K,'Value',sliceIndexK);
set(jSlider_T,'HighValue',Tf(2));
set(jSlider_T,'LowValue',Tf(1));

setThresholdFunc([],[],{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K});

end

%% Plot slices
function plotSlice(~,~,inputCell)

hf=inputCell{1};
jSlider=inputCell{2};
dirOpt=inputCell{3};
sliceIndex = get(jSlider,'Value');
hf.UserData.sliceIndices(dirOpt)=sliceIndex;
sliceIndices=hf.UserData.sliceIndices;

M=hf.UserData.M;
logicThreshold=hf.UserData.logicThreshold;
v=hf.UserData.v;
patchType=hf.UserData.patchTypes{dirOpt};

logicPatch=false(size(M));
switch dirOpt
    case 1
        logicPatch(sliceIndex,:,:)=1;
    case 2
        logicPatch(:,sliceIndex,:)=1;
    case 3
        logicPatch(:,:,sliceIndex)=1;
end
logicPatch=logicPatch & logicThreshold;

[F,V,C]=ind2patch(logicPatch,M,patchType);
[V(:,1),V(:,2),V(:,3)]=im2cart(V(:,2),V(:,1),V(:,3),v);

if isfield(hf.UserData,'hp')
    try
        delete(hf.UserData.hp(dirOpt));
    catch
    end
end

hf.UserData.hp(dirOpt)= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','none','FaceAlpha',hf.UserData.faceAlpha);
if dirOpt==3
    if hf.UserData.showAll<0
        plotContourSet(hf);
    end
end
navString=['I: ',num2str(sliceIndices(1)),', J:  ',num2str(sliceIndices(2)),', K: ',num2str(sliceIndices(3))];
hf.UserData.colorBarhandle.Label.String=navString;
hf.Name=[hf.UserData.Name,' ',navString];
drawnow;

end

%% Change threshold

function setThresholdFunc(~,~,inputCell)

hf=inputCell{1};
jSlider_T=inputCell{2};
jSlider_I=inputCell{3};
jSlider_J=inputCell{4};
jSlider_K=inputCell{5};

Tf(1) = get(jSlider_T,'LowValue');
Tf(2) = get(jSlider_T,'HighValue');

M=hf.UserData.M;
W=max(M(:))-min(M(:));

T_low=min(M(:))+(W*Tf(1)/100);
T_high=min(M(:))+(W*Tf(2)/100);
logicThreshold=(M>=T_low & M<=T_high);
hf.UserData.logicThreshold=logicThreshold;

plotSlice([],[],{hf,jSlider_I,1});
plotSlice([],[],{hf,jSlider_J,2});
plotSlice([],[],{hf,jSlider_K,3});

end

%% Scroll bar resizing

function setScrollSizeFunc(~,~,inputCell)
hf=inputCell{1};
w=inputCell{2};

for q=3:6
    jSlider=inputCell{q};
    javacomponent(jSlider,[w*(q-3),0,w,round(hf.Position(4))]);
end

% hColorBar=hf.UserData.colorBarhandle;
% hColorBar.Units='pixels';
% hColorBar.AxisLocation='in';
% hColorBar.Position(1)=hf.Position(3)-50;
% hColorBar.Position(2)=50;
% hColorBar.Position(4)=hf.Position(4)-100;

hColorBar=hf.UserData.colorBarhandle;
hColorBar.Units='pixels';
hColorBar.AxisLocation='in';
hColorBar.Position(1)=4*60;
hColorBar.Position(2)=10;
hColorBar.Position(3)=hf.Position(3)-(4*60)-10;
hColorBar.Position(4)=25;

% hPopUp1=hf.UserData.hPopUp1;
% hPopUp1.Position=[hf.Position(3)-100 hf.Position(4)-50 100 50];

hAxis=hf.UserData.hAxis;
axis_units=hAxis.Units;
hAxis.Units='pixels';
hAxis.Position(2)=50;
hAxis.Position(4)=hf.Position(4)-100;
hAxis.Units=axis_units;

end

%% Help

function helpFunc(~,~)

msgText={'INPUT OPTIONS:',...
    '------------------------------------------------------------------------',...
    'Up              -> Increase contour level',...
    'Down            -> Decrease contour level',...
    'Left click      -> Select closest contour to keep',...
    'Right click     -> Delete closest contour',...
    'C               -> Define crop window for contours',...
    'D               -> Manually draw points as a contour, left click to draw, right to finish',...
    'V               -> Add selected contours from the previous slice for use in current slice',...
    'O               -> Reorder contour points. Warning this is a slow process work in ordered e.g. clockwise fashion to avoid ordering',...
    'R               -> Remove all selected/drawn contours',...
    'S               -> Create split contour, i.e. keep current contour segment but add aditional for this slice',...
    'Space bar       -> Done with current slice, proceed to next',...
    'H               -> Display help',...
    'M               -> Alter colorbar limits',...
    'L               -> Jump to contour level',...
    'J               -> Change contour increment size (increment when up key is pressed)',...
    'Esc               -> Force finish and exit',...
    };
helpdlg(msgText,'Help information');

end

%% Save

function saveFunc(~,~,inputCell)
hf=inputCell{1}; %Figure handle

%Dialog for save path and save name
prompt = {'Save path (leave empty to browse to desired folder instead):','Save name (no extension):'};
dlg_title = 'Saving';
defaultOptions = {hf.UserData.savePath,hf.UserData.saveName};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);

if ~isempty(Q)
    %Deal with empty path
    if isempty(Q{1})
        Q{1}=uigetdir(hf.UserData.savePath,'Select save path'); %Browse
        if Q{1}==0
            return;
        end
    end
    
    %Create save folder if it does not exist already
    if ~exist(Q{1},'dir')
        mkdir(Q{1});
    end
    
    %Override UserData struct
    hf.UserData.savePath=Q{1};
    hf.UserData.saveName=Q{2};
    
    %Save UserData
    
    saveFields={'sliceIndices',...
        'contourSlice',...
        'fontColor',...
        'faceAlpha',...
        'lineColors',...
        'sketchContour',...
        'ContourSet',...
        'colormapSketch',...
        'colormapSet',...
        'csapsSmoothPar',...
        'showAll',...
        'logicThreshold'};
    
    for q=1:1:numel(saveFields)
        fieldNameCurrent=saveFields{q};
        saveStruct.(fieldNameCurrent)=hf.UserData.(fieldNameCurrent);
    end
    
    save(fullfile(hf.UserData.savePath,[hf.UserData.saveName,'.mat']),'saveStruct');
end

end

%% Load

function loadFunc(~,~,inputCell)
hf=inputCell{1}; %Figure handle

%Dialog for save path and save name
prompt = {'File (leave empty to browse to desired file instead):'};
dlg_title = 'Loading';
defaultOptions = {fullfile(hf.UserData.savePath,hf.UserData.saveName)};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);

if ~isempty(Q)
    %Deal with empty path
    if isempty(Q{1})
        [loadName, loadPath, ~] = uigetfile({'*.mat','MAT-files (*.mat)'},'Pick a file','MultiSelect', 'off');
        [~,loadName,~] = fileparts(loadName); %Remove extension
    else
        [loadPath,loadName,~] = fileparts(Q{1});
    end
    
    %Override UserData struct
    hf.UserData.savePath=loadPath;
    hf.UserData.saveName=loadName;
    
    %Import content
    load(fullfile(loadPath,[loadName,'.mat']))
    try
        fieldSet = fieldnames(saveStruct); % Cell containing all structure field names
        for q=1:1:numel(fieldSet)
            fieldNameCurrent=fieldSet{q};
            hf.UserData.(fieldNameCurrent)=saveStruct.(fieldNameCurrent);
        end
    catch
        hf.UserData.ContourSet=Vcs; 
    end
    
    %Update plots
    plotSketchContour(hf);
    plotContourSet(hf);
    
    %Change sliders
    set(hf.UserData.sliderHandles{1},'Value',hf.UserData.sliceIndices(1));
    set(hf.UserData.sliderHandles{2},'Value',hf.UserData.sliceIndices(2));
    set(hf.UserData.sliderHandles{3},'Value',hf.UserData.sliceIndices(3));
    
    %Update slices
    plotSlice([],[],{hf,hf.UserData.sliderHandles{1},1});
    plotSlice([],[],{hf,hf.UserData.sliderHandles{2},2});
    plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
    
end

end

%% Home

function homeFunc(~,~,inputCell)
hf=inputCell{1}; %Figure handle
set(hf.UserData.sliderHandles{3},'Value',hf.UserData.contourSlice);
plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
view(0,90);
hAxis=hf.UserData.hAxis;
axis_units=hAxis.Units;
hAxis.Units='pixels';
hAxis.Position(2)=50;
hAxis.Position(4)=hf.Position(4)-100;
hAxis.Units=axis_units;
end

%% colormap
% 
% function funcColormap(source,~)
% selectIndex = source.Value;
% cMapStrings = source.String;
% cMapStringNew = cMapStrings{selectIndex};
% colormap(cMapStringNew);
% end

%% colorbar

function colorbarFunc(~,~)
prompt = {'Minimum:','Maximum:'};
dlg_title = 'Set colorbar limits';

defaultOptions = {'0','1'};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    minC=str2double(Q{1});
    maxC=str2double(Q{2});
    caxis([minC maxC]);
end

end

%% ghost

function ghostFunc(~,~,inputCell)

hf=inputCell{1}; %Figure handle

prompt = {'Enter alpha level:'};
dlg_title = 'Set transparency';

defaultOptions = {'0.8'};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    hf.UserData.faceAlpha=str2double(Q{1});
    
    %Update slices
    plotSlice([],[],{hf,hf.UserData.sliderHandles{1},1});
    plotSlice([],[],{hf,hf.UserData.sliderHandles{2},2});
    plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
end

end

%% showHide

function showHideFunc(~,~,inputCell)
hf=inputCell{1}; %Figure handle
hf.UserData.showAll=hf.UserData.showAll*-1; %Change state
%Update plot
plotContourSet(hf);
end

%% draw

function drawFunc(~,~,inputCell)

hf=inputCell{1}; %Figure handle
view(0,90); 
[qSlice]=updateSliceIndex(hf);

v=hf.UserData.v;
zs=(qSlice-0.5).*v(3); %z-level for contourslice


hpd=[];
drawDone=0;
Vd=[];
[mousePointerType]=specialPointerShape('pen');
while drawDone==0
    [Xd,Yd,bd]=qginput(1,mousePointerType);
    switch bd
        case 1
            Vd=[Vd; Xd Yd zs];
            delete(hpd);
            hpd=plotV(Vd,'b.-','MarkerSize',20,'lineWidth',2);
        case 3
            drawDone=1;
    end
end
delete(hpd);

if ~isempty(Vd)
    
    D=pathLength(Vd);
    n=2*round(max(D(:))./min(v(1:2)));
    [Vd] = evenlySampleCurve(Vd,n,'pchip',0);
    
    %Get contours
    C=hf.UserData.sketchContour;
    
    %Add drawing to contours
    C{end+1}=Vd;
    
    %Override sketch contour
    hf.UserData.sketchContour=C;
    
    %Update plot
    plotSketchContour(hf);    
    
end

%Turn off button
set(hf.UserData.ButtonHandles.Draw,'State','Off');

setDefaultPointer; %Set default pointer

end

%% cut

function cutFunc(~,~,inputCell)

hf=inputCell{1}; %Figure handle

set(hf,'WindowButtonDownFcn','','WindowButtonUpFcn','');

v=hf.UserData.v;
siz=size(hf.UserData.M);
qSlice=hf.UserData.contourSlice;

zs=(qSlice-0.5).*v(3); %z-level for contourslice

homeFunc([],[],inputCell);

%Distance thresholds for curve spacing, higher spacings are considered gaps for instance
distThresh=2*max(v(1:2));

%Using rbbox to draw a cropping window
[mousePointerType]=specialPointerShape('cut');
set(gcf,'Pointer','custom','PointerShapeCData',mousePointerType.PointerShapeCData,'PointerShapeHotSpot',mousePointerType.PointerShapeHotSpot); %User specified mousePointerType

k = waitforbuttonpress;
if k==0 %If mouse click
    
    point1 = get(gca,'CurrentPoint');    % button down detected
    rbbox;                 % return figure units
    drawnow; %pause(1e-6); % Required to avoid point1=point2 errors
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);            % extract x and y
    point2 = point2(1,1:2);
    p1 = [min(point1(1,1),point2(1,1)) max(point1(1,2),point2(1,2))];
    p2 = [max(point1(1,1),point2(1,1)) min(point1(1,2),point2(1,2))];
    [ic1,jc1,~]=cart2im(p1(1),p1(2),0,v);
    [ic2,jc2,~]=cart2im(p2(1),p2(2),0,v);
    hcc1=plot3(p1(1),p1(2),zs,'r+','MarkerSize',20);
    hcc2=plot3(p2(1),p2(2),zs,'r+','MarkerSize',20);
    drawnow;
    
    mc=false(siz(1:2));
    mc(round(ic2):round(ic1),round(jc1):round(jc2))=1;
    cropInd=find(mc);
    
    %Get contours
    C=hf.UserData.sketchContour;
    
    %Remove contour elements within crop window and regroup
    num_C=numel(C);
    L_delete=false(1,num_C);
    for ig=1:1:num_C
        V=C{ig};
        
        [icg,jcg,~]=cart2im(V(:,1),V(:,2),zeros(size(V(:,1))),v);
        icg=round(icg);
        jcg=round(jcg);
        icg(icg<1)=1;
        jcg(jcg<1)=1;
        icg(icg>siz(1))=siz(1);
        jcg(jcg>siz(2))=siz(2);
        
        vertexInd= sub2ind(siz(1:2),icg,jcg);
        L_crop=~ismember(vertexInd,cropInd);
        
        if any(~L_crop(:))
            if all(L_crop(:))
                L_delete(ig)=1; %Whole contour will be removed
            else %Contour is cropped
                %Split curve in groups
                groupIndices=cumsum(diff([0; L_crop]).*L_crop).*L_crop;
                groupIndices(groupIndices==0)=NaN;
                
                D_startEnd=sqrt(sum((V(1,:)-V(end,:)).^2));
                if all(~ismember([1 size(V,1)],cropInd)) && D_startEnd<distThresh% If the start and end are not cropped reunite start/end groups
                    groupIndices(groupIndices==1)=max(groupIndices(:));
                end
                
                groupFirstIter=1;
                for q_group=nanmin(groupIndices(:)):1:nanmax(groupIndices(:))
                    
                    Vg=V(groupIndices==q_group,:); %Current curve
                    
                    %Reorder so within curve gaps are removed
                    Dx=diff(Vg(:,1)); Dy=diff(Vg(:,2));
                    D=sqrt(Dx.^2+Dy.^2);
                    LD=D>distThresh;
                    if any(LD) %reorder
                        indGap=find(LD,1)+1;
                        Vg=[Vg(indGap:end,:); Vg(1:indGap-1,:)];
                    end
                    
                    if groupFirstIter==1
                        C{ig}=Vg;
                        groupFirstIter=0;
                    else
                        C{end+1}=Vg;
                    end
                end
            end
        end
    end
    
    indDelete=find(L_delete);
    L_delete=false(1,numel(C));
    L_delete(indDelete)=1;
    
    C=C(~L_delete);
    
    %Get new group sizes
    %Get group sizes
    group_sizes = cell2mat(cellfun(@(x) size(x,1), C,'UniformOutput',0)');
    
    %Remove contours of lenght 1
    C=C(group_sizes>1);
    group_sizes=group_sizes(group_sizes>1);
    
    %Sort C array according to group size
    [group_sizes,ind_sort]=sort(group_sizes);
    C=C(ind_sort);
    
    delete(hcc1); delete(hcc2);
    
    %Override sketch contour
    hf.UserData.sketchContour=C;
    
    %Update plot
    plotSketchContour(hf);
    
    
end

setDefaultPointer; %Set default pointer

%Turn off button
set(hf.UserData.ButtonHandles.Cut,'State','Off');

set(hf,'KeyPressFcn', {@figKeyPressFunc,{hf}},'WindowButtonDownFcn', {@figMouseDown,{hf}},'WindowButtonUpFcn', {@mouseup,hf});

end

%% sample

function sampleFunc(~,~,inputCell)

hf=inputCell{1}; %Figure handle

view(0,90);

v=hf.UserData.v;
M=hf.UserData.M;
X=hf.UserData.X;
Y=hf.UserData.Y;
Z=hf.UserData.Z;

qSlice=hf.UserData.sliceIndices(3); %Get current slice
hf.UserData.contourSlice=qSlice; %Set current slice as the contour slice

zs=(qSlice-0.5).*v(3); %z-level for contourslice

while 1
    [xc,yc,b]=qginput(1);
    if ~isempty(b)
        switch b
            case 1
                vClick=[xc yc zs];
                [i,j,~]=cart2im(xc,yc,zs,v);
                i=round(i);
                j=round(j);
                k=qSlice;
                Tc=M(i,j,k);
                
                %Compute contour
                [C]=gcontour(X(:,:,qSlice),Y(:,:,qSlice),M(:,:,qSlice),Tc,min(v)/4,'pchip');
                
                %Get the closest contour to the click
                [indMin]=findNearestContour(C,vClick);
                                
                V=C{indMin};
                
                V(:,3)=zs*ones(size(V,1),1);
                
                hf.UserData.sketchContour={V};
                
                %Update plot
                plotSketchContour(hf);
            otherwise
                break
        end
    else
        break
    end
end
set(hf.UserData.ButtonHandles.Sample,'State','Off');

end

%% Figure key press

function figKeyPressFunc(src,eventData,inputCell)

hf=inputCell{1}; %Figure handle

% step = 1;
% if ismember('shift', eventData.Modifier)
%     step = -step; %Make negative while shift is down
% end
%
% if ismember('control', eventData.Modifier)
%     step = step * 4; %Increase speed
% end

% Key input options
switch eventData.Key
    case 'leftarrow' % Pan left
        a=hf.UserData.faceAlpha;
        a=(a-0.1);
        if a<0
            a=0;
        end
        hf.UserData.faceAlpha=a;
        
        %Update slices
        plotSlice([],[],{hf,hf.UserData.sliderHandles{1},1});
        plotSlice([],[],{hf,hf.UserData.sliderHandles{2},2});
        plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
    case 'rightarrow' % Pan right
        a=hf.UserData.faceAlpha;
        a=(a+0.1);
        if a>1
            a=1;
        end
        hf.UserData.faceAlpha=a;
        
        %Update slices
        plotSlice([],[],{hf,hf.UserData.sliderHandles{1},1});
        plotSlice([],[],{hf,hf.UserData.sliderHandles{2},2});
        plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
    case 'downarrow' % Pan down
        qSlice=hf.UserData.sliceIndices(3);
        if qSlice>1
            hf.UserData.sliceIndices(3)=hf.UserData.sliceIndices(3)-1;
            set(hf.UserData.sliderHandles{3},'Value',hf.UserData.sliceIndices(3));
            plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
            if hf.UserData.showAll<0
                plotContourSet(hf);
            end
        end
    case 'uparrow' % Pan up
        qSlice=hf.UserData.sliceIndices(3);
        if qSlice<size(hf.UserData.M,3)
            hf.UserData.sliceIndices(3)=hf.UserData.sliceIndices(3)+1;
            set(hf.UserData.sliderHandles{3},'Value',hf.UserData.sliceIndices(3));
            plotSlice([],[],{hf,hf.UserData.sliderHandles{3},3});
             if hf.UserData.showAll<0
                plotContourSet(hf);
            end
        end
    case 'v' % Activate vcw
        set(hf.UserData.cFigure.Handles.vcw,'State','On');
    case 's' % Sample
        set(hf.UserData.ButtonHandles.Sample,'State','On');
        sampleFunc([],[],inputCell);
    case 'c' % Cut
        set(hf.UserData.ButtonHandles.Cut,'State','On');
        cutFunc([],[],inputCell);
    case 'home' % Home
        homeFunc([],[],inputCell);
    case 'd'
        set(hf.UserData.ButtonHandles.Draw,'State','On');
        drawFunc([],[],inputCell);
    case 'a'
        set(hf.UserData.ButtonHandles.Select,'State','On');
        selectFunc([],[],inputCell);
    case 'q'
        set(hf.UserData.ButtonHandles.hSmooth,'State','On');
        smoothFunc([],[],inputCell);
    case 'delete'
        set(hf.UserData.ButtonHandles.Delete,'State','On');
        deleteFunc([],[],inputCell);
    case 'equal'
        set(hf.UserData.ButtonHandles.hGrow,'State','On');
        inputCell{2}=1;
        growShrinkFunc([],[],inputCell);
    case 'hyphen'
        set(hf.UserData.ButtonHandles.hShrink,'State','On');
        inputCell{2}=-1;
        growShrinkFunc([],[],inputCell);
end

end

%% Figure click

function figMouseDown(src, eventData,inputCell)

hf=inputCell{1};

funcs = {'pan','rot','zoomz','zoomz'};

% Get the button pressed
% cax = overobj2('axes');

cax = get(hf, 'CurrentAxes');
if isempty(cax)
    return;
end

checkAxisLimits(hf);
colorbarLocSet(hf,'manual');

switch get(hf, 'SelectionType')
    case 'extend' % Middle button
        mouseDownFunc = ['vcw_',funcs{2}];
    case 'alt' % Right hand button
        mouseDownFunc = ['vcw_',funcs{3}];
    case 'open' % Double click
        caxUserDataStruct=get(cax,'UserData');
        camview(cax,caxUserDataStruct.defaultView);
        return;
    otherwise
        mouseDownFunc = ['vcw_',funcs{1}];
end

% Set the cursor
switch mouseDownFunc
    case {'vcw_zoom', 'vcw_zoomz'}
        shape=[ 2   2   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   1   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   2 NaN NaN NaN   2   2   2   2  ;
            2   1   2   1   1   2   1   1   1   2 NaN   2   1   2   1   2  ;
            2   1   2   1   2 NaN   2   1   1   1   2   1   1   2   1   2  ;
            2   2   2   2 NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   1   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   2   2  ];
    case 'vcw_pan'
        shape=[ NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
            NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
            2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
            2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
            NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
            NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ];
    case {'vcw_rotz', 'vcw_rot'}
        % Rotate
        shape=[ NaN NaN NaN   2   2   2   2   2 NaN   2   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN   1   1   1   1   1   2   1   1   2 NaN NaN NaN NaN ;
            NaN NaN NaN   2   1   1   1   1   2   1   1   1   2 NaN NaN NaN ;
            NaN NaN   2   1   1   1   1   1   2   2   1   1   1   2 NaN NaN ;
            NaN   2   1   1   1   2   1   1   2 NaN NaN   2   1   1   2 NaN ;
            NaN   2   1   1   2 NaN   2   1   2 NaN NaN   2   1   1   2 NaN ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            NaN   2   1   1   2 NaN NaN   2   1   2 NaN   2   1   1   2 NaN ;
            NaN   2   1   1   2 NaN NaN   2   1   1   2   1   1   1   2 NaN ;
            NaN NaN   2   1   1   1   2   2   1   1   1   1   1   2 NaN NaN ;
            NaN NaN NaN   2   1   1   1   2   1   1   1   1   2 NaN NaN NaN ;
            NaN NaN NaN NaN   2   1   1   2   1   1   1   1   1 NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   2 NaN   2   2   2   2   2 NaN NaN NaN ];
    otherwise
        return
end
mouseDownFunc=str2func(mouseDownFunc); 

% Record where the pointer is
global VCW_POS
VCW_POS = get(0, 'PointerLocation');

% Set the cursor and callback
set(hf, 'Pointer', 'custom', 'pointershapecdata', shape, 'WindowButtonMotionFcn', {mouseDownFunc, cax});

end


%%
function mouseup(src, eventData,hf)
% Clear the cursor and callback
set(hf, 'WindowButtonMotionFcn', '', 'Pointer', 'arrow');
end


%%
function d = check_vals(s, d)
% Check the inputs to the manipulation methods are valid
global VCW_POS
if ~isempty(s)
    % Return the mouse pointers displacement
    new_pt = get(0, 'PointerLocation');
    d = VCW_POS - new_pt;
    VCW_POS = new_pt;
end
end

%% Figure manipulation functions
function vcw_rot(s, d, cax, hf)
d = check_vals(s, d);
try
    % Rotate XY
    camorbit(cax, d(1), d(2), 'camera', [0 0 1]);
catch
    % Error, so release mouse down
    mouseup(hf);
end
end

function vcw_rotz(s, d, cax, hf)
d = check_vals(s, d);
try
    % Rotate Z
    camroll(cax, d(2));
catch
    % Error, so release mouse down
    mouseup(hf);
end
end

function vcw_zoom(s, d, cax, hf)
d = check_vals(s, d);
% Zoom
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2));
try
    camzoom(cax, d);
catch
    % Error, so release mouse down
    mouseup(hf);
end
end

function vcw_zoomz(s, d, cax, hf)
d = check_vals(s, d);
% Zoom by moving towards the camera
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2)) - 1;
try
    camdolly(cax, 0, 0, d, 'fixtarget', 'camera');
catch
    % Error, so release mouse down
    mouseup(hf);
end
end

function vcw_pan(s, d, cax, hf)
d = check_vals(s, d);
try
    % Pan
    camdolly(cax, d(1), d(2), 0, 'movetarget', 'pixels');
catch
    % Error, so release mouse down
    mouseup(hf);
end
end

function colorbarLocSet(hf,locOpt)
H=findobj(hf,'Type','colorbar'); %Handle set
for q=1:1:numel(H)
    if isa(locOpt,'cell')
        set(H(q),'Location',locOpt{q});
    else
        set(H(q),'Location',locOpt);
    end
end
end


function checkAxisLimits(hf)

h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
if ~isempty(h)
    for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
        axis(h);
        
        xLim=get(h,'xlim');
        yLim=get(h,'ylim');
        zLim=get(h,'zlim');
        
        wx=abs(diff(xlim));
        wy=abs(diff(ylim));
        wz=abs(diff(zlim));
        
        w_max=max([wx wy wz]);
        min_w=1e-3;
        if w_max<min_w
            w_max=min_w;
        end
        
        w_min=w_max/10;
        if w_min<min_w
            w_min=min_w;
        end
        w_add=[-w_min w_min]/2;       
        
        if wx<w_min
            set(h,'xlim',xLim+w_add);
        end
        
        if wy<w_min
            set(h,'ylim',yLim+w_add);
        end
        
        if wz<w_min
            set(h,'zlim',zLim+w_add);
        end              
    end
    drawnow; 
end

end


%% Select contour

function selectFunc(~,~,inputCell)

homeFunc([],[],inputCell);

hf=inputCell{1}; %Figure handle
v=hf.UserData.v;

%Get contours
C=hf.UserData.sketchContour;

if ~isempty(C)
    
    homeFunc([],[],inputCell);
    [mousePointerType]=specialPointerShape('smallHand');
    qSlice=hf.UserData.contourSlice; %Get contour slice    
    zs=(qSlice-0.5).*v(3); %z-level for contourslice
    
    while 1
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    Vcs=hf.UserData.ContourSet;
                    
                    vClick=[xc yc zs];
                    
                    %Find closest contour
                    [indMin]=findNearestContour(C,vClick);
                    V_pick=C{indMin};
                    
                    if isempty(Vcs{qSlice}{1})
                        Vcs{qSlice}{1}=V_pick; %Add as first new
                    else
                        %Ask user what to do with this contour
                        A = questdlg('What would you like to do with this contour?','Adding contour','New','Merge','Replace','Replace');
                        
                        if ~isempty(A)
                            switch A
                                case 'New'
                                    Vcs{qSlice}{end+1}=V_pick;
                                case 'Replace'
                                    [xc,yc,b]=qginput(1,mousePointerType);
                                    if ~isempty(b)
                                        switch b
                                            case 1
                                                vClick=[xc yc zs];
                                                [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                                                Vcs{qSlice}{indMin_Vcs}=V_pick; %Replace
                                            otherwise
                                                break
                                        end
                                    end
                                case 'Merge' %Add selected contour to contour set
                                    
                                    [xc,yc,b]=qginput(1,mousePointerType);
                                    if ~isempty(b)
                                        switch b
                                            case 1
                                                vClick=[xc yc zs];
                                                [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                                                
                                                %Reorder selection before adding to maintain overall
                                                %contour point order.
                                                
                                                %See if start or end is closest to current contour end
                                                V_now=Vcs{qSlice}{indMin_Vcs}; %The contour so far
                                                
                                                V_now_start_end=V_now([1 size(V_now,1)],:);
                                                V_add_start_end=V_pick([1 size(V_pick,1)],:);
                                                try
                                                    D=dist(V_now_start_end,V_add_start_end');
                                                catch
                                                    D=distND(V_now_start_end,V_add_start_end);
                                                end
                                                [~,indMin_lin]=min(D(:));
                                                
                                                switch indMin_lin
                                                    case 1 %1 1 -> starts are close, add in front and flip
                                                        V_pick=flipud(V_pick);
                                                        Vd=[V_pick; V_now];
                                                    case 2 %2 1 -> start is close to end, just add behind
                                                        Vd=[V_now; V_pick];
                                                    case 3 %1 2 -> end is close to start, just add in fron
                                                        Vd=[V_pick; V_now];
                                                    case 4 %2 2 -> ends are clsoe, add behind and flip
                                                        V_pick=flipud(V_pick);
                                                        Vd=[V_now; V_pick];
                                                end
                                                
                                                D=pathLength(Vd);
                                                n=2*round(max(D(:))./min(v(1:2)));
                                                [Vd] = evenlySampleCurve(Vd,n,'pchip',0);
                                                
                                                Vcs{qSlice}{indMin_Vcs}=Vd;
                                                
                                            otherwise
                                                break
                                        end
                                    end
                            end
                        else
                            break
                        end
                    end
                    %Override contour set
                    hf.UserData.ContourSet=Vcs;
                    
                    %Update plot
                    plotContourSet(hf);
                    
                    %Delete from draft contour list
                    L_keep=true(1,numel(C)); L_keep(indMin)=0;
                    C=C(L_keep);
                    
                    %Override sketch contour
                    C=fixEmptyContours(C);
                    hf.UserData.sketchContour=C;
                    
                    
                    %Update plot
                    plotSketchContour(hf);
                    
                otherwise
                    break
            end
        else
            break
        end
        if isempty(C)
            break
        end
    end
else
    warndlg('No contours assigned to current slice');
end

set(hf.UserData.ButtonHandles.Select,'State','Off');
setDefaultPointer;

end

%% Close
function closeFunc(~,~,inputCell)

hf=inputCell{1};
view(0,90);
[qSlice]=updateSliceIndex(hf);

v=hf.UserData.v;
zs=(qSlice-0.5).*v(3); %z-level for contourslice
Vcs=hf.UserData.ContourSet;

if ~isempty(Vcs{qSlice}{1})
    [mousePointerType]=specialPointerShape('smallHand');
    while 1        
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    vClick=[xc yc zs];                    
                    [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                    if ~isempty(Vcs{qSlice}{indMin_Vcs})
                        Vd=Vcs{qSlice}{indMin_Vcs};
                        D=pathLength([Vd;Vd(1,:)]);
                        n=2*round(max(D(:))./min(v(1:2)));
                        [Vd] = evenlySampleCurve(Vd,n,'pchip',1);
                        Vcs{qSlice}{indMin_Vcs}=Vd;
                    end
                    hf.UserData.ContourSet{qSlice}=Vcs{qSlice};
                    plotContourSet(hf);  
                otherwise
                    break
            end
        end
    end
end

end

%% Smooth

function smoothFunc(~,~,inputCell)

hf=inputCell{1};
view(0,90);
[qSlice]=updateSliceIndex(hf);

v=hf.UserData.v;
zs=(qSlice-0.5).*v(3); %z-level for contourslice
Vcs=hf.UserData.ContourSet;

hf.UserData.sketchContour=fixEmptyContours(hf.UserData.sketchContour);
numSketchedContours=numel(hf.UserData.sketchContour);
if numSketchedContours==1
    if isempty(hf.UserData.sketchContour{1})
        numSketchedContours=0; 
    end
end

if ~isempty(Vcs{qSlice}{1})
    [mousePointerType]=specialPointerShape('smallHand');
    while 1        
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    vClick=[xc yc zs];                    
                    [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                    
                    Vd=Vcs{qSlice}{indMin_Vcs}; %Get closest contour
                    
                    hf.UserData.sketchContour{numSketchedContours+1}=Vd; %Add unsmoothened contour as sketched contour
                                         
                    D=pathLength([Vd;Vd(1,:)]);
                    n=2*round(max(D(:))./min(v(1:2)));
                    p=hf.UserData.csapsSmoothPar;
                    [Vd] = evenlySampleCurve(Vd,n,p,1);
                    Vcs{qSlice}{indMin_Vcs}=Vd;
                    hf.UserData.ContourSet{qSlice}=Vcs{qSlice};                                        

                    plotContourSet(hf);
                    plotSketchContour(hf);                 
                    
                case 112 % p                    
                    prompt = {'Enter smoothening parameter:'};
                    dlg_title = 'Adjust smoothening parameter';
                    defaultOptions = {num2str(hf.UserData.csapsSmoothPar)};
                    s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
                    Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);                    
                    if ~isempty(Q)
                        hf.UserData.csapsSmoothPar=str2double(Q{1});
                    end
                otherwise
                    break
            end
        end
    end
end

set(hf.UserData.ButtonHandles.hSmooth,'State','Off');
setDefaultPointer;

end

%% Ellipse fit

function ellipseFunc(~,~,inputCell)

hf=inputCell{1};
view(0,90);
[qSlice]=updateSliceIndex(hf);

v=hf.UserData.v;
zs=(qSlice-0.5).*v(3); %z-level for contourslice
Vcs=hf.UserData.ContourSet;

if ~isempty(Vcs{qSlice}{1})
    [mousePointerType]=specialPointerShape('smallHand');
    while 1        
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    vClick=[xc yc zs];                    
                    [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                    
                    Vd=Vcs{qSlice}{indMin_Vcs}; %Get closest contour
                                                           
                    [A] = ellipseFit(Vd,2,[]);             

                    t=linspace(0,2*pi,size(Vd,1));
                    t=t(1:end-1);                    
                    [Vd]=ellipseCoord(A,t);
                    D=pathLength([Vd;Vd(1,:)]);
                    n=2*round(max(D(:))./min(v(1:2)));
                    [Vd] = evenlySampleCurve(Vd,n,'pchip',1);
                    Vd(:,3)=zs;
                    Vcs{qSlice}{indMin_Vcs}=Vd;
                    
                    hf.UserData.ContourSet{qSlice}=Vcs{qSlice};                                        

                    plotContourSet(hf);
                    
                otherwise
                    break
            end
        end
    end
end

set(hf.UserData.ButtonHandles.hEllipse,'State','Off');
setDefaultPointer;

end

%% Delete sketched contour

function deleteFunc(~,~,inputCell)

hf=inputCell{1};
v=hf.UserData.v;
%Get contours
C=hf.UserData.sketchContour;
homeFunc([],[],inputCell);
[mousePointerType]=specialPointerShape('x');
qSlice=hf.UserData.contourSlice; %Get contour slice
zs=(qSlice-0.5).*v(3); %z-level for contourslice

while 1
    if ~isempty(C)
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    vClick=[xc yc zs];
                    
                    %Find closest contour
                    [indMin]=findNearestContour(C,vClick);
                    L_keep=true(1,numel(C)); L_keep(indMin)=0;
                    C=C(L_keep);
                    
                    %Override sketch contour
                    C=fixEmptyContours(C);
                    hf.UserData.sketchContour=C;
                    
                    %Update plot
                    plotSketchContour(hf);
                otherwise
                    break
            end
        else
            break
        end
    else
        break
    end
end
setDefaultPointer;
set(hf.UserData.ButtonHandles.Delete,'State','Off');

end

%% reset

function resetFunc(~,~,inputCell)

hf=inputCell{1}; %Figure handle

[qSlice]=updateSliceIndex(hf);

hf.UserData.ContourSet{qSlice}={[]};
plotContourSet(hf);

% hf.UserData.sketchContour={};
% plotSketchContour(hf);

end

%% Convert contour

function convertFunc(~,~,inputCell)

homeFunc([],[],inputCell);

hf=inputCell{1};

qSlice=hf.UserData.sliceIndices(3); %Get current slice
hf.UserData.contourSlice=qSlice;

v=hf.UserData.v;
zs=(qSlice-0.5).*v(3); %z-level for contourslice

%Get current sketch contours
C=hf.UserData.sketchContour;

%Get current contour set
Vcs=hf.UserData.ContourSet;

[mousePointerType]=specialPointerShape('smallHand');
while 1
    if ~isempty(Vcs{qSlice}{1})
    [xc,yc,b]=qginput(1,mousePointerType);
    if ~isempty(b)
        switch b
            case 1                                
                vClick=[xc yc zs];
                [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                V_pick=Vcs{qSlice}{indMin_Vcs}; %Contour to downgrade
                Vcs{qSlice}{indMin_Vcs}=[];
                Vcs{qSlice}=fixEmptyContours(Vcs{qSlice});
                hf.UserData.ContourSet{qSlice}=Vcs{qSlice};
                
%                 if isempty(C{1})
%                     C{1}=V_pick;
%                 else
                    C{end+1}=V_pick;
%                 end
                C=fixEmptyContours(C);
                
                hf.UserData.sketchContour=C;
                
                %Update plots
                plotSketchContour(hf);
                plotContourSet(hf);                
            otherwise
                break
        end
    end
    else        
        break
    end
end

end

%% Grow
function growShrinkFunc(~,~,inputCell)

growDir=inputCell{2};
hf=inputCell{1};
view(0,90);
[qSlice]=updateSliceIndex(hf);

view(0,90);

v=hf.UserData.v;
M=hf.UserData.M;
X=hf.UserData.X;
Y=hf.UserData.Y;
Z=hf.UserData.Z;

zs=(qSlice-0.5).*v(3); %z-level for contourslice
Vcs=hf.UserData.ContourSet;


mask_height=1;
mask_shape=[mask_height:-1:1 -1:-1:-mask_height]';
MASK_I=[mask_shape; zeros(size(mask_shape));]';
MASK_J=[zeros(size(mask_shape)); mask_shape;]';
MASK_K=[zeros(size(mask_shape)); zeros(size(mask_shape));]';


if ~isempty(Vcs{qSlice}{1})
    [mousePointerType]=specialPointerShape('smallHand');
    while 1        
        [xc,yc,b]=qginput(1,mousePointerType);
        if ~isempty(b)
            switch b
                case 1
                    vClick=[xc yc zs];                    
                    [indMin_Vcs]=findNearestContour(Vcs{qSlice},vClick);
                    
                    Vd=Vcs{qSlice}{indMin_Vcs}; %Get closest contour
                     
                    L = inpolygon(X(:,:,qSlice),Y(:,:,qSlice),Vd(:,1),Vd(:,2));                    
                    x=X(:,:,qSlice);
                    y=Y(:,:,qSlice);
                    [D,~]=minDist([x(:) y(:)],Vd(:,[1 2]));                    
                    D(L)=-D(L);
                    m=reshape(D,[size(M,1) size(M,2)]);
                    
                    if growDir==1
                        Tc=min(v/4);
                    elseif growDir==-1
                        Tc=-min(v/4);
                    end
                    
                    %Compute contour
                    [c,group_sizes]=gcontour(X(:,:,qSlice),Y(:,:,qSlice),m,Tc,min(v)/4,'pchip');
                    
                    [~,maxInd]=max(group_sizes);
                    Vd=c{maxInd};
                    Vd(:,3)=zs*ones(size(Vd,1),1);
                    
                    Vcs{qSlice}{indMin_Vcs}=Vd;
                    
                    hf.UserData.ContourSet{qSlice}=Vcs{qSlice};                                        

                    plotContourSet(hf);
                    
                otherwise
                    break
            end
        end
    end
end

if growDir==1
    set(hf.UserData.ButtonHandles.hGrow,'State','Off');
elseif growDir==-1
    set(hf.UserData.ButtonHandles.hShrink,'State','Off');
end

setDefaultPointer;

end

%% Plot sketch slice contours

function plotSketchContour(hf)

%Remove current sketch contour if it exists
if isfield(hf.UserData,'sketchContourHandle')
    if ~isempty(hf.UserData.sketchContourHandle)
        delete(hf.UserData.sketchContourHandle);
        hf.UserData.sketchContourHandle=[];
    end
else
    hf.UserData.sketchContourHandle=[];
end

%Plot new
C=hf.UserData.sketchContour;
plotColors=resampleColormap(hf.UserData.colormapSketch,numel(C));%autumn(numel(C));
for q=1:1:numel(C)    
    V=C{q};
    if ~isempty(V)
        h=plotV(V,'r--','lineWidth',3);
        set(h,'Color',plotColors(q,:));
        hf.UserData.sketchContourHandle=[hf.UserData.sketchContourHandle h];
    end
end

end

%% Plot contour set

function plotContourSet(hf)

%Remove current sketch contour if it exists
if isfield(hf.UserData,'sketchContourHandle')
    if ~isempty(hf.UserData.contourSetHandle)
        delete(hf.UserData.contourSetHandle);
        hf.UserData.contourSetHandle=[];
    end
else
    hf.UserData.contourSetHandle=[];
end

%Plot new
Vcs=hf.UserData.ContourSet;

if hf.UserData.showAll<0
    sliceRange=updateSliceIndex(hf);
%     sliceRange=hf.UserData.sliceIndices(3);
else    
    sliceRange=1:1:numel(Vcs);        
end

for qSlice=sliceRange
    for qc=1:1:numel(Vcs{qSlice})
        plotColors=resampleColormap(hf.UserData.colormapSet,numel(Vcs{qSlice}));
        if ~isempty(Vcs{qSlice}{qc})
            V=Vcs{qSlice}{qc};
            h=plotV(V,'g-','lineWidth',3); set(h,'Color',plotColors(qc,:));
            hf.UserData.contourSetHandle=[hf.UserData.contourSetHandle h];
        end
    end
end
drawnow;

end

%% Setting default pointer

function setDefaultPointer
set(gcf,'Pointer','arrow');%'watch');
end

%% Find nearest sketch contour

function [indMin]=findNearestContour(C,vClick)
D=NaN(1,numel(C));
for qC=1:1:numel(C)
    V=C{qC};    
    if ~isempty(V)    
        D(qC)=min(sqrt(sum((V(:,[1 2])-vClick(ones(size(V,1),1),[1 2])).^2,2)));
    end
end
[~,indMin] = min(D);
end

%% Fix empty contours

function C=fixEmptyContours(C)

numEntries = cell2mat(cellfun(@(x) numel(x), C,'UniformOutput',0)');
logicKeep=numEntries>0;
if all(~logicKeep)
    C={[]};
else
    C=C(logicKeep);
end

end

%% Update slice index

function [qSlice]=updateSliceIndex(hf)

qSliceContour=hf.UserData.contourSlice; %Get current slice index
qSlice=hf.UserData.sliceIndices(3); %Get new slice index

if qSlice~=qSliceContour
    hf.UserData.contourSlice=qSlice;
    hf.UserData.sketchContour={}; 
    plotSketchContour(hf);    
end

end
