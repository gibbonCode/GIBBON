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

		UserData = get(hf,'UserData');

        %Add view control widget
		if is_octave
			% To do: vcw handling
		else
			hp=vcw(hf);
			UserData.cFigure.Handles.vcw=hp;
		end


        %Add export figure widget
		if is_octave
			% To do: efw handling
		else
			efw(hf);
		end

        if ~isfield(UserData,'anim8')
            error('Loaded figure lacks an anim8 structure as UserData');
        else
            animStruct=UserData.anim8.animStruct; %Get anim8 structure
        end

        try %Try to remove a slider if present
            hSlider=UserData.anim8.sliderHandles{1};
			delete(hSlider)
        catch %Nothing to delete
        end

        try %Try to remove text label above slider if present
            delete(hf.UserData.anim8.ButtonHandles.hTextTime)
        catch %Nothing to delete
        end
		set(hf,'UserData',UserData);
    case 2 %Create new
        hf=varargin{1}; %Figure handle
        animStruct=varargin{2}; %The anim8 structure
		UserData = get(hf,'UserData');
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

animTime = animStruct.Time(:);
numSliderSteps = numel(animTime);

sliceIndexI = numSliderSteps; %Initial index at end
sliderStep = [1/(numSliderSteps-1) 1/(numSliderSteps-1)]; %Slider step sizes

%Initialize slider
p = get(hf,'Position');
hSlider= uicontrol(hf,'Style','slider','Position',[0,0,round(p(3)),scrollBarWidth]);
set(hSlider,'Value',sliceIndexI,'Min',1,'Max',numSliderSteps,'SliderStep',sliderStep);
set(hSlider,'Callback',{@updateViewFunc,hf});
% hSlider.KeyPressFcn={@updateViewFunc,hf};

if is_octave
	% Continuous monitoring of the slider is perhaps not implemented yet in Octave
	% Basic discrete (once mouse is released) listener
	addlistener(hSlider,'Value',@(hObject, event) updateViewFunc(hObject, event,hf));
else
	addlistener(hSlider,'ContinuousValueChange',@(hObject, event) updateViewFunc(hObject, event,hf));
	addlistener(hSlider,'Value','PostSet',@(hObject, event) updateViewFunc(hObject, event,hf));
end

%


%% Set resize function
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
    hb = findall(hf,'Type','');

    %Remove them if they are not the default bar
    if ~isempty(hb)
        for q=1:1:numel(hb)
			tagString = get(hb(q),'Tag');
			% The following checks for the default MATLAB and Octave toolbar tag strings respectively
            if ~any(strcmp(tagString,{'FigureToolBar','__default_toolbar__'}));
                delete(hb(q));
            end
        end
    end
end

%Create new anim8 toolbar
hb = uitoolbar(hf);
set(hb,'Tag',toolbarTag);

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
S(isnan(S))=1; % Dealing with transparency bug in Octave
iconHelp=S;

% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Help','CData',iconHelp,'Tag','help_button','ClickedCallback',@helpFunc,'Separator','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Play button
% This button will switch appearance from a play icon to a stop icon

%get icon 1
D=importdata(fullfile(iconPath,'play.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
S(isnan(S))=1; % Dealing with transparency bug in Octave
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
S(isnan(S))=1; % Dealing with transparency bug in Octave
iconStop=S;

% Create a uitoggletool in the toolbar
hPlay=uitoggletool(hb,'TooltipString','Play','CData',iconPlay,'Tag','play_button','ClickedCallback',{@playFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time button

%get icon 3
D=importdata(fullfile(iconPath,'time.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
S(isnan(S))=1; % Dealing with transparency bug in Octave
iconTime=S;

% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Time','CData',iconTime,'Tag','time_button','ClickedCallback',{@timeFunc,{hf}}); %,'Separator','on'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cycle button

%get icon 4
D=importdata(fullfile(iconPath,'cycle.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
S(isnan(S))=1; % Dealing with transparency bug in Octave
iconCycle=S;

% Create a uitoggletool in the toolbar
hCycle=uitoggletool(hb,'TooltipString','Cycle forward-backward','CData',iconCycle,'Tag','cycle_button');

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
S(isnan(S))=1; % Dealing with transparency bug in Octave
iconSave=S;

% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Save .gif animation','CData',iconSave,'Tag','saveAnimation_button','ClickedCallback',{@saveAnimationFunc,{hf}}); %,'Separator','on'

%% Text fields

% title(sprintf('%6.16e ',T),'FontSize',hf.UserData.anim8.fontSize);
% title(sprintf('%f',T),'FontSize',hf.UserData.anim8.fontSize);
figPos = get(hf,'Position');
hTextTime = uicontrol(hf,'Style','text',...
    'String',[' Time: ',sprintf('%f',animTime(1))],...
    'Position',[0 scrollBarWidth round(figPos(3)) scrollBarWidth],...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','Left',...
    'fontSize',fontSize);

%% Modify figure UserData
UserData.anim8.sliderHandles={hSlider};
% hf.UserData.anim8.ButtonHandles.Sample=hSample;
UserData.anim8.animStruct=animStruct;
UserData.anim8.fontSize=fontSize;
% animHandles=animStruct.Handles; %Handles of objects to animate
% animProps=animStruct.Props; %Properties of objects to animate
% animSets=animStruct.Set; %Property values for to set in order to animate
% animTime=animStruct.Time;
UserData.anim8.playDir=1;
UserData.anim8.ButtonHandles.Play=hPlay;
UserData.anim8.ButtonHandles.hCycle=hCycle;
UserData.anim8.ButtonHandles.hTextTime=hTextTime;
UserData.anim8.pauseTime=1/numel(animStruct.Time);%mean(diff(animStruct.Time));
UserData.anim8.shiftMag=1;%ceil(numel(animStruct.Time)/20);
UserData.anim8.icons.play=iconPlay;
UserData.anim8.icons.stop=iconStop;

% Store current settings
UserData.anim8.WindowButtonDownFcn = get(hf,'WindowButtonDownFcn');
UserData.anim8.WindowButtonUpFcn = get(hf,'WindowButtonUpFcn');
UserData.anim8.KeyPressFcn = get(hf,'KeyPressFcn');
UserData.anim8.WindowScrollWheelFcn = get(hf,'WindowScrollWheelFcn');
UserData.anim8.BusyAction=get(hf,'BusyAction');

%Export figure widget settings
UserData.efw.defaultPath=fullfile(cd,'efw');
UserData.efw.imName=['figure',num2str(get(hf,'Number'))];
UserData.efw.imExt='png';
UserData.efw.imRes='50';
UserData.efw.exportFigOpt='-nocrop';
UserData.efw.exportGifOpt='1';

set(hf,'UserData',UserData);

%% Set initial slider position and force update of graphics
% Set at arbitrary state that is different from the desired initial
if sliceIndexI<numSliderSteps
	set(hSlider,'Value',sliceIndexI+1); % Add one if not at max
else
	set(hSlider,'Value',sliceIndexI-1); % Subtract one if at max
end
set(hSlider,'Value',sliceIndexI); % Now set at desired initial
drawnow;

end

%% Scroll bar resizing

function figResize(~,~,inputCell)
hf = inputCell{1};
scrollBarWidth = inputCell{2};
UserData = get(hf,'UserData');

p = round(get(hf,'Position'));
hSlider = inputCell{3};

set(hSlider,'Position',[0,0,p(3), scrollBarWidth]);
set(UserData.anim8.ButtonHandles.hTextTime,'Position',[0 scrollBarWidth p(3) scrollBarWidth]);

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
hf = inputCell{1};
UserData = get(hf,'UserData');
shiftMag=UserData.anim8.shiftMag;

set(UserData.anim8.ButtonHandles.Play,'CData',UserData.anim8.icons.stop,'TooltipString','Stop');

while strcmp(get(UserData.anim8.ButtonHandles.Play,'State'),'on')
    tic
    if strcmp(get(UserData.anim8.ButtonHandles.hCycle,'State'),'on')
        hSlider=UserData.anim8.sliderHandles{1};

        sliderValue=get(hSlider,'Value');
        sliderValueNew=sliderValue+(shiftMag*UserData.anim8.playDir);

        sliderMax=get(hSlider,'Max');
        sliderMin=get(hSlider,'Min');

        if sliderValueNew<sliderMin
            UserData.anim8.playDir=UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*UserData.anim8.playDir);
			set(hf,'UserData',UserData);
        elseif sliderValueNew>sliderMax
            UserData.anim8.playDir=UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*UserData.anim8.playDir);
			set(hf,'UserData',UserData);
        end
        set(hSlider,'Value',sliderValueNew);
    else
        shiftSlider(UserData.anim8.sliderHandles{1},shiftMag);
    end
    updateViewFunc([],[],hf);
    drawnow;

    t=toc;
    tPause=UserData.anim8.pauseTime-t;
    if tPause>0
        pause(tPause);
    end
end

set(UserData.anim8.ButtonHandles.Play,'CData',UserData.anim8.icons.play,'TooltipString','Play');

end

%% Time
function timeFunc(~,~,inputCell)
hf = inputCell{1};
UserData = get(hf,'UserData');
prompt = {'Enter pause time:','Enter step size:'};
dlg_title = 'Time stepping settings';
defaultOptions = {num2str(UserData.anim8.pauseTime),num2str(UserData.anim8.shiftMag)};
s=40+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s; 1 s;],defaultOptions);
if ~isempty(Q)
    if ~isempty(Q{1})
        UserData.anim8.pauseTime=str2double(Q{1});
    end
    if ~isempty(Q{2})
        UserData.anim8.shiftMag=str2double(Q{2});
    end
	set(hf,'UserData',UserData);
end

end

%% Save Animation

function saveAnimationFunc(~,~,inputCell)

hf = inputCell{1};
UserData = get(hf,'UserData');

defStruct = UserData.efw; %Default option structure
exportGifAnim8(hf,defStruct,1); %Export gif data

end

%% Figure key press

function figKeyPressFunc(~,eventData,inputCell)

hf = inputCell{1}; %Figure handle
UserData = get(hf,'UserData');

% Key input options
switch eventData.Key
    case {'leftarrow','downarrow'}
        shiftSlider(UserData.anim8.sliderHandles{1},-1*UserData.anim8.shiftMag);
        updateViewFunc([],[],hf);
    case {'rightarrow','uparrow'}
        shiftSlider(UserData.anim8.sliderHandles{1},1*UserData.anim8.shiftMag);
        updateViewFunc([],[],hf);
    case 'home' % Go to the start
        hSlider=UserData.anim8.sliderHandles{1};
        sliderVal=get(hSlider,'Min');
        set(hSlider,'Value',sliderVal);
        updateViewFunc([],[],hf);
    case 'end' % Go to the end
        hSlider=UserData.anim8.sliderHandles{1};
        sliderVal=get(hSlider,'Max');
        set(hSlider,'Value',sliderVal);
        updateViewFunc([],[],hf);
    case 'v' % Activate vcw
        set(hf.UserData.cFigure.Handles.vcw,'State','On');
    case'space'
        if strcmp(get(UserData.anim8.ButtonHandles.Play,'State'),'off')
            set(UserData.anim8.ButtonHandles.Play,'State','on');
            playFunc([],[],{hf});
        else
            set(UserData.anim8.ButtonHandles.Play,'State','off');
        end
end

end

%% updateViewFunc

function updateViewFunc(~,~,hf)

UserData = get(hf,'UserData');

% Get slider handle and set value
hSlider = UserData.anim8.sliderHandles{1}; %Get the handle from the structure
sliderValue = round(get(hSlider,'Value')); %Round just in case, since we need whole numbers
set(hSlider,'Value',sliderValue); %Set slider value

% Get time, handles, properties, and property values to set
T = UserData.anim8.animStruct.Time(sliderValue); % Current slider "time" stamp
H = UserData.anim8.animStruct.Handles{sliderValue};%e.g. [hp,hp]; %Handles of objects to animate
P = UserData.anim8.animStruct.Props{sliderValue};% e.g. {'FaceColor','Vertices'}; %Properties of objects to animate
S = UserData.anim8.animStruct.Set{sliderValue};%e.g.{c(q,:),V}; %Property values for to set in order to animate

% Update the title string
set(UserData.anim8.ButtonHandles.hTextTime,'String',[' Time: ',sprintf('%f',T)],'FontSize',UserData.anim8.fontSize);

% Update the properties for this slider position
for q=1:1:numel(H)
    h = H(q); % Current graphics handle
    p = P{q}; % Current graphics property
    s = S{q}; % Current property setting
    set(h,p,s);% Setting the property
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
