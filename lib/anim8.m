function hf=anim8(varargin)

%% Parse input

switch nargin
    case 2
        hf=varargin{1};
        animStruct=varargin{2};        
end

%% Visualisation settings

% fontColor='w';
fontSize=15;
% cMap=gjet(250);

scrollBarWidth=50; 

%% Defining slider 

animTime=animStruct.Time(:);

minT=1;
maxT=numel(animTime);
w=maxT-minT;
sliceIndexI=1; %Initial index
% nTickMajor=round(w/20);
tickSizeMajor_I=round(w/20);

%% Initialize display

figure(hf);

%Initialize slider
jSlider = javax.swing.JSlider(minT,maxT);
javacomponent(jSlider,[0,0,round(hf.Position(3)),scrollBarWidth]);
set(jSlider, 'MajorTickSpacing',tickSizeMajor_I, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@updateViewFunc,{hf,jSlider}},'Orientation',jSlider.HORIZONTAL);

%% Set resize function

set(hf,'ResizeFcn',{@setScrollSizeFunc,{hf,scrollBarWidth,jSlider}});

%% Initialize figure callbacks

set(hf,'KeyPressFcn', {@figKeyPressFunc,{hf}},'WindowButtonDownFcn', {@figMouseDown,{hf}},'WindowButtonUpFcn', {@mouseup,hf});

%% Initialise buttons

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
iconPath=fullfile(toolboxPath,'icons');

% hb = findall(hf,'Type','uitoolbar');
hb = uitoolbar(hf);

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

%get icon
D=importdata(fullfile(iconPath,'play.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uitoggletool in the toolbar
hPlay=uitoggletool(hb,'TooltipString','Play forward','CData',S,'Tag','play_button','ClickedCallback',{@playFunc,{hf}});

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

%% Text fields

% title(sprintf('%6.16e ',T),'FontSize',hf.UserData.fontSize);
% title(sprintf('%f',T),'FontSize',hf.UserData.fontSize);

hTextTime = uicontrol(hf,'Style','text',...
                'String',[' Time: ',sprintf('%f',animTime(1))],...
                'Position',[0 scrollBarWidth round(hf.Position(3)) round(scrollBarWidth/2)],'BackgroundColor',[1 1 1],'HorizontalAlignment','Left'); %hf.Color
            
%% Set figure UserData

hf.UserData.sliderHandles={jSlider};
% hf.UserData.ButtonHandles.Sample=hSample;
hf.UserData.animStruct=animStruct;
hf.UserData.fontSize=fontSize;
% animHandles=animStruct.Handles; %Handles of objects to animate
% animProps=animStruct.Props; %Properties of objects to animate
% animSets=animStruct.Set; %Property values for to set in order to animate
% animTime=animStruct.Time;
hf.UserData.playDir=1;
hf.UserData.ButtonHandles.Play=hPlay;
hf.UserData.ButtonHandles.hCycle=hCycle;
hf.UserData.ButtonHandles.hTextTime=hTextTime;
hf.UserData.pauseTime=2/numel(animStruct.Time);
hf.UserData.shiftMag=1; 

% Store current settings
hf.UserData.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.BusyAction=hf.BusyAction;

%% Initialize slider locations
set(jSlider,'Value',sliceIndexI);

end

%% Scroll bar resizing

function setScrollSizeFunc(~,~,inputCell)
hf=inputCell{1};
scrollBarWidth=inputCell{2};

jSlider=inputCell{3};
javacomponent(jSlider,[0,0,round(hf.Position(3)), scrollBarWidth]);

set(hf.UserData.ButtonHandles.hTextTime,'Position',[0 scrollBarWidth round(hf.Position(3)) round(scrollBarWidth/2)]);
end

%% Help

function helpFunc(~,~)

msgText={'INPUT OPTIONS:',...
    '------------------------------------------------------------------------',...
    'Up              -> Increase contour level',...    
    };
helpdlg(msgText,'Help information');

end

%% play

function playFunc(~,~,inputCell)
hf=inputCell{1};
shiftMag=hf.UserData.shiftMag; 

while strcmp(get(hf.UserData.ButtonHandles.Play,'State'),'on')
    tic
    if strcmp(get(hf.UserData.ButtonHandles.hCycle,'State'),'on')
        
        jSlider=hf.UserData.sliderHandles{1};
        
        sliderValue=get(jSlider,'Value');
        sliderValueNew=sliderValue+(shiftMag*hf.UserData.playDir);
        
        sliderMax=get(jSlider,'Maximum');
        sliderMin=get(jSlider,'Minimum');
        
        if sliderValueNew<sliderMin
            hf.UserData.playDir=hf.UserData.playDir*-1;            
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.playDir);
        elseif sliderValueNew>sliderMax
            hf.UserData.playDir=hf.UserData.playDir*-1;            
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.playDir);
        end
                
        set(jSlider,'Value',sliderValueNew);
        
    else
        shiftSlider(hf.UserData.sliderHandles{1},shiftMag);
    end
    
    
    drawnow;
    
    t=toc;
    if (hf.UserData.pauseTime-t)>0        
        pause(hf.UserData.pauseTime-t);
    end    
end

end

%% Time
function timeFunc(~,~,inputCell)
hf=inputCell{1};

prompt = {'Enter pause time:','Enter step size'};
dlg_title = 'Pause time (s)';
defaultOptions = {num2str(hf.UserData.pauseTime),num2str(hf.UserData.shiftMag)};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    if ~isempty(Q{1})
        hf.UserData.pauseTime=str2double(Q{1});
    end
    if ~isempty(Q{2})
        hf.UserData.shiftMag=str2double(Q{2});
    end
end

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
    case {'leftarrow','downarrow'}
         shiftSlider(hf.UserData.sliderHandles{1},-1*hf.UserData.shiftMag);
    case {'rightarrow','uparrow'}
        shiftSlider(hf.UserData.sliderHandles{1},1*hf.UserData.shiftMag);
    case 'v' % Activate vcw
        set(hf.UserData.cFigure.Handles.vcw,'State','On');
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

%% Setting default pointer

function setDefaultPointer
set(gcf,'Pointer','arrow');%'watch');
end

%% updateViewFunc

function updateViewFunc(src,eventData,inputCell)

hf=inputCell{1};
jSlider=inputCell{2};
animStruct=hf.UserData.animStruct;

sliderValue=get(jSlider,'Value');

T=animStruct.Time(sliderValue); 
H=animStruct.Handles{sliderValue};%e.g. [hp,hp]; %Handles of objects to animate
P=animStruct.Props{sliderValue};% e.g. {'FaceColor','Vertices'}; %Properties of objects to animate
S=animStruct.Set{sliderValue};%e.g.{c(q,:),V}; %Property values for to set in order to animate

set(hf.UserData.ButtonHandles.hTextTime,'String',[' Time: ',sprintf('%f',T)],'FontSize',hf.UserData.fontSize);

for q=1:1:numel(H)    
    h=H(q); % Current graphics handle
    p=P{q}; % Current graphics property
    s=S{q}; %Current property setting
    h.(p)=s;
end

end

%% Shift slider

function shiftSlider(jSlider,shiftMag)

sliderValue=get(jSlider,'Value');
sliderValueNew=sliderValue+shiftMag;

sliderMax=get(jSlider,'Maximum'); 
sliderMin=get(jSlider,'Minimum'); 

if sliderValueNew<sliderMin
    set(jSlider,'Value',sliderMax);
elseif sliderValueNew>sliderMax
    set(jSlider,'Value',1);
else
    set(jSlider,'Value',sliderValueNew);
end

end

