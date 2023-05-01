function hf=anim8(varargin)

% function hf=anim8(hf,animStruct)
% ------------------------------------------------------------------------
% 
% Change log: 
% 2019/05/20 Added loading of saved anim8 figure by using the path to the
% figure as sole input (or no input which triggers uigetfile)
% 2019/08/09 Changed to use uicontrol slider rather than java slider due
% to future removal of javacomponent
% 2019/10/13 Fixed bug in relation to handling a single time step (no
% annimation), fixed by copying state twice and animating clones states. 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case {0,1} %Load previous 
        if nargin==1 %If path to figure is specified 
            figureFullPath=varargin{1};
        else %If no path is specified
            [figureFileName,figurePath]=uigetfile;
            figureFullPath=fullfile(figurePath,figureFileName);
        end
        hf=open(figureFullPath); %Open figure and keep handle
          
        %Add view control widget
        hp=vcw(hf);
        hf.UserData.cFigure.Handles.vcw=hp;
        
        %Add export figure widget
        efw(hf);
        
        if ~isfield(hf.UserData,'anim8')
            error('Loaded figure lacks an anim8 structure as UserData');
        else
            animStruct=hf.UserData.anim8.animStruct; %Get anim8 structure
        end
        
        try %Try to remove a slider if present
            hSlider=hf.UserData.anim8.sliderHandles{1};
            if isa(hSlider,'matlab.ui.control.UIControl')
                delete(hSlider)
            else %Old jSlider type 
                try %javacomponent is scheduled to be removed but try it here for old slider bar type
                    warning off %Disable all warnings
                    [~,hc]=javacomponent(hSlider);
                    delete(hc);
                    warning on %Enable warnings again
                catch
                    %No fix if the above throws an error
                end
            end
        catch %Nothing to delete
        end
        
        try %Try to remove text label above slider if present
            delete(hf.UserData.anim8.ButtonHandles.hTextTime)
        catch %Nothing to delete
        end        
        
    case 2 %Create new
        hf=varargin{1}; %Figure handle
        animStruct=varargin{2}; %The anim8 structure
end

if ~isfield(animStruct,'Time')
    animStruct.Time=linspace(0,numel(animStruct.Handles),numel(animStruct.Handles));
end

%% Visualisation settings

% fontColor='w';
fontSize=15;
% cMap=gjet(250);
scrollBarWidth=30;

%% Defining slider
figure(hf); drawnow;

if numel(animStruct.Time)==1
    warning('Only 1 step for animation');
    animStruct.Time(end+1)=animStruct.Time(end);
    animStruct.Handles{end+1}=animStruct.Handles{end}; %Handles of objects to animate
    animStruct.Props{end+1}=animStruct.Props{end}; %Properties of objects to animate
    animStruct.Set{end+1}=animStruct.Set{end};
end

animTime=animStruct.Time(:);
sliceIndexI=numel(animTime); %Initial index at end
sliderStep=[1/(numel(animTime)-1) 1/(numel(animTime)-1)]; %Slider step sizes

%Initialize slider
hSlider= uicontrol(hf,'Style','slider','Position',[0,0,round(hf.Position(3)),scrollBarWidth]);
set(hSlider,'Value',sliceIndexI,'Min',1,'Max',numel(animTime),'SliderStep',sliderStep);
hSlider.Callback={@updateViewFunc,hf};
% hSlider.KeyPressFcn={@updateViewFunc,hf};
addlistener(hSlider,'ContinuousValueChange',@(hObject, event) updateViewFunc(hObject, event,hf));
addlistener(hSlider,'Value','PostSet',@(hObject, event) updateViewFunc(hObject, event,hf));

%% Set resize function

% set(hf,'ResizeFcn',{@figResize,{hf,scrollBarWidth,hSlider}});
% set(hf,'ResizeFcn',@(h,e)figResize(h,e,{hf,scrollBarWidth,hSlider}));

hFunc=get(hf,'ResizeFcn');

if iscell(hFunc)
    warning('anim8 changed the ResizeFcn function for this figure. Specify your ResizeFcn in the form @(h,e)figResize(h,e,c) to avoid this behavior.');    
    set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,scrollBarWidth,hSlider}));
else
    if isempty(hFunc)
        set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,scrollBarWidth,hSlider}));
    else        
        set(hf,'ResizeFcn',@(a,b)(cellfun(@(x)feval(x,a,b),{hFunc,@(a,b)figResize(a,b,{hf,scrollBarWidth,hSlider})})));
    end
end


%% Initialize figure callbacks

set(hf,'KeyPressFcn', {@figKeyPressFunc,{hf}});%'WindowButtonDownFcn', {@figMouseDown,{hf}},'WindowButtonUpFcn', {@mouseup,hf});

%% Initialise tool bar and buttons

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
iconPath=fullfile(toolboxPath,'icons');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Toolbar
toolbarTag='anim8_toolbar';

%Remove existing anim8 toolbar
hb = findall(hf,'Tag',toolbarTag); % hb = findall(hf,'Type','uitoolbar');
if ~isempty(hb)     
    delete(hb)
else    
    %Find all toolbars
    hb = findall(hf,'Type','uitoolbar');
    
    %Remove them if they are not the default bar
    if ~isempty(hb)
        for q=1:1:numel(hb)
            if ~strcmp(hb(q).Tag,'FigureToolBar')
                delete(hb(q));
            end
        end
    end    
end

%Create new anim8 toolbar
hb = uitoolbar(hf);
hb.Tag=toolbarTag;

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
%Play button

%get icon 1
D=importdata(fullfile(iconPath,'play.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
iconPlay=S;

%get icon 2
D=importdata(fullfile(iconPath,'stop.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
iconStop=S;

% Create a uitoggletool in the toolbar
hPlay=uitoggletool(hb,'TooltipString','Play','CData',iconPlay,'Tag','play_button','ClickedCallback',{@playFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time button

%get icon
D=importdata(fullfile(iconPath,'time.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Time','CData',S,'Tag','time_button','ClickedCallback',{@timeFunc,{hf}}); %,'Separator','on'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cycle button

%get icon
D=importdata(fullfile(iconPath,'cycle.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uitoggletool in the toolbar
hCycle=uitoggletool(hb,'TooltipString','Cycle forward-backward','CData',S,'Tag','cycle_button');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save animation button

%get icon
D=importdata(fullfile(iconPath,'saveAnimation.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Save .gif animation','CData',S,'Tag','saveAnimation_button','ClickedCallback',{@saveAnimationFunc,{hf}}); %,'Separator','on'

%% Text fields

% title(sprintf('%6.16e ',T),'FontSize',hf.UserData.anim8.fontSize);
% title(sprintf('%f',T),'FontSize',hf.UserData.anim8.fontSize);

hTextTime = uicontrol(hf,'Style','text',...
    'String',[' Time: ',sprintf('%f',animTime(1))],...
    'Position',[0 scrollBarWidth round(hf.Position(3)) scrollBarWidth],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left',...
    'fontSize',fontSize);

%% Set figure UserData

hf.UserData.anim8.sliderHandles={hSlider};
% hf.UserData.anim8.ButtonHandles.Sample=hSample;
hf.UserData.anim8.animStruct=animStruct;
hf.UserData.anim8.fontSize=fontSize;
% animHandles=animStruct.Handles; %Handles of objects to animate
% animProps=animStruct.Props; %Properties of objects to animate
% animSets=animStruct.Set; %Property values for to set in order to animate
% animTime=animStruct.Time;
hf.UserData.anim8.playDir=1;
hf.UserData.anim8.ButtonHandles.Play=hPlay;
hf.UserData.anim8.ButtonHandles.hCycle=hCycle;
hf.UserData.anim8.ButtonHandles.hTextTime=hTextTime;
hf.UserData.anim8.pauseTime=1/numel(animStruct.Time);%mean(diff(animStruct.Time));
hf.UserData.anim8.shiftMag=1;%ceil(numel(animStruct.Time)/20);
hf.UserData.anim8.icons.play=iconPlay;
hf.UserData.anim8.icons.stop=iconStop;

% Store current settings
hf.UserData.anim8.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.anim8.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.anim8.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.anim8.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.anim8.BusyAction=hf.BusyAction;

%Export figure widget settings
hf.UserData.efw.defaultPath=fullfile(cd,'efw');
hf.UserData.efw.imName=['figure',num2str(get(hf,'Number'))];
hf.UserData.efw.imExt='png';
hf.UserData.efw.imRes='50';
hf.UserData.efw.exportFigOpt='-nocrop';
hf.UserData.efw.exportGifOpt='1';

%% Initialize slider locations
set(hSlider,'Value',round(sliceIndexI/2));
set(hSlider,'Value',sliceIndexI);
drawnow;

end

%% Scroll bar resizing

function figResize(~,~,inputCell)
hf=inputCell{1};
scrollBarWidth=inputCell{2};

hSlider=inputCell{3};
if isa(hSlider,'matlab.ui.control.UIControl') %Check if its a MATLAB style slider
    set(hSlider,'Position',[0,0,round(hf.Position(3)), scrollBarWidth]);
    set(hf.UserData.anim8.ButtonHandles.hTextTime,'Position',[0 scrollBarWidth round(hf.Position(3)) scrollBarWidth]);    
else %Old JAVA jSlider
%     try
% %         delete(hSlider)
%     catch
%     end
%     try %javacomponent is scheduled to be removed but try it here for old slider bar type
%         jc=javacomponent(hSlider,[0,0,round(hf.Position(3)), scrollBarWidth]);
%         jc.delete
%     catch
%         %No fix if the above throws an error
%     end
end

end

%% Help

function helpFunc(~,~)

msgText={'* Play button or press the space bar -> Start animation',...
    '* Stop button or press the space bar -> Stop animation',...
    '* Time button -> Change time stepping settings',...
    '* Cycle button -> Toggle forward / forward-backward mode',...
    '* Save button -> Export snap shorts (e.g. jpg) and animated .gif files',...
    '* Press the v key to activate the view control widget (vcw), note that vcw changes key press functions until vcw is deactivated',...
    };
helpdlg(msgText,'Help information');

end

%% play

function playFunc(~,~,inputCell)
hf=inputCell{1};
shiftMag=hf.UserData.anim8.shiftMag;

set(hf.UserData.anim8.ButtonHandles.Play,'CData',hf.UserData.anim8.icons.stop,'TooltipString','Stop');

while strcmp(get(hf.UserData.anim8.ButtonHandles.Play,'State'),'on')
    tic
    if strcmp(get(hf.UserData.anim8.ButtonHandles.hCycle,'State'),'on')        
        hSlider=hf.UserData.anim8.sliderHandles{1};
        
        sliderValue=get(hSlider,'Value');
        sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        
        sliderMax=get(hSlider,'Max');
        sliderMin=get(hSlider,'Min');
        
        if sliderValueNew<sliderMin
            hf.UserData.anim8.playDir=hf.UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        elseif sliderValueNew>sliderMax
            hf.UserData.anim8.playDir=hf.UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        end        
        set(hSlider,'Value',sliderValueNew);        
    else
        shiftSlider(hf.UserData.anim8.sliderHandles{1},shiftMag);
    end
    updateViewFunc([],[],hf);
    drawnow;
    
    t=toc;
    tPause=hf.UserData.anim8.pauseTime-t; 
    if tPause>0
        pause(tPause);
    end
end

set(hf.UserData.anim8.ButtonHandles.Play,'CData',hf.UserData.anim8.icons.play,'TooltipString','Play');

end

%% Time
function timeFunc(~,~,inputCell)
hf=inputCell{1};

prompt = {'Enter pause time:','Enter step size:'};
dlg_title = 'Time stepping settings';
defaultOptions = {num2str(hf.UserData.anim8.pauseTime),num2str(hf.UserData.anim8.shiftMag)};
s=40+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    if ~isempty(Q{1})
        hf.UserData.anim8.pauseTime=str2double(Q{1});
    end
    if ~isempty(Q{2})
        hf.UserData.anim8.shiftMag=str2double(Q{2});
    end
end

end

%% Save Animation

function saveAnimationFunc(~,~,inputCell)

hf=inputCell{1};
defStruct=hf.UserData.efw; %Default option structure
exportGifAnim8(hf,defStruct,1); %Export gif data

end

%% Figure key press

function figKeyPressFunc(~,eventData,inputCell)

hf=inputCell{1}; %Figure handle
    
% Key input options
switch eventData.Key
    case {'leftarrow','downarrow'}
        shiftSlider(hf.UserData.anim8.sliderHandles{1},-1*hf.UserData.anim8.shiftMag);
        updateViewFunc([],[],hf);
    case {'rightarrow','uparrow'}
        shiftSlider(hf.UserData.anim8.sliderHandles{1},1*hf.UserData.anim8.shiftMag);
        updateViewFunc([],[],hf);
    case 'home' % Go to the start
        hSlider=hf.UserData.anim8.sliderHandles{1};
        sliderVal=get(hSlider,'Min');
        set(hSlider,'Value',sliderVal);        
        updateViewFunc([],[],hf);
    case 'end' % Go to the end
        hSlider=hf.UserData.anim8.sliderHandles{1};
        sliderVal=get(hSlider,'Max');
        set(hSlider,'Value',sliderVal);
        updateViewFunc([],[],hf);
    case 'v' % Activate vcw
        set(hf.UserData.cFigure.Handles.vcw,'State','On');
    case'space'
        if strcmp(get(hf.UserData.anim8.ButtonHandles.Play,'State'),'off')
            set(hf.UserData.anim8.ButtonHandles.Play,'State','on');
            playFunc([],[],{hf});
        else
            set(hf.UserData.anim8.ButtonHandles.Play,'State','off');
        end
end

end

%% updateViewFunc

function updateViewFunc(~,~,hf)

animStruct=hf.UserData.anim8.animStruct;
hSlider=hf.UserData.anim8.sliderHandles{1};
sliderValue=round(get(hSlider,'Value'));
set(hSlider,'Value',sliderValue);

T=animStruct.Time(sliderValue);
H=animStruct.Handles{sliderValue};%e.g. [hp,hp]; %Handles of objects to animate
P=animStruct.Props{sliderValue};% e.g. {'FaceColor','Vertices'}; %Properties of objects to animate
S=animStruct.Set{sliderValue};%e.g.{c(q,:),V}; %Property values for to set in order to animate

set(hf.UserData.anim8.ButtonHandles.hTextTime,'String',[' Time: ',sprintf('%f',T)],'FontSize',hf.UserData.anim8.fontSize);

for q=1:1:numel(H)
    h=H(q); % Current graphics handle
    p=P{q}; % Current graphics property
    s=S{q}; % Current property setting    
    h.(p)=s;% Setting the property
end

end

%% Shift slider

function shiftSlider(hSlider,shiftMag)

sliderValue=get(hSlider,'Value');
sliderValueNew=sliderValue+shiftMag;

sliderMax=get(hSlider,'Max');
sliderMin=get(hSlider,'Min');

if sliderValueNew<sliderMin
    set(hSlider,'Value',sliderMax);
elseif sliderValueNew>sliderMax
    set(hSlider,'Value',1);
else
    set(hSlider,'Value',sliderValueNew);
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
