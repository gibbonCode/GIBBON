function [varargout]=vcw(varargin)

% function vcw(hf,buttonOpt)
% ------------------------------------------------------------------------
% vcw, View Control Widget, Allows users to manipulate a view in 3D using
% key and button presses
%
%   vcw(hf,buttonOpt)
%
% Allows the user to rotate, pan and zoom a figure using key presses and
% mouse gestures. 
%
% Input:
%   hf - Handle of the figure to be manipulated (default: gcf).
%   buttonOpt - 4x1 cell array indicating the function to associate with
%             each mouse button (left to right) and the scroll action.
%             Functions can be any of:
%                'rot' - Rotate about x and y axes of viewer's coordinate
%                        frame
%                'rotz' - Rotate about z axis of viewer's coordinate frame
%                'zoom' - Zoom (change camera view angle)
%                'zoomz' - Move along z axis of viewer's coordinate frame
%                'pan' - Pan
%                '' - Don't use that button
%             Default: {'pan','rot','zoomz','zoomz'};).
%
% This code was inspired by the fcw function by Oliver Woodford (which was
% based on Torsten Vogel's view3d function, which was in turn inspired by
% rotate3d from The MathWorks, Inc.). Although vcw is inspired by and
% similar to fcw function by Oliver Woodford, it offers many added features
% and bug-fixes.  
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/04/15 Copied from fcw and renamed to vcw due to planned revision
% 2015/04/15 Added: 1) handing of colorbars (bug in fcw when view(2) is
% used combined with panning which induced zooming and panning), 2) overobj
% axes selection so that the current axes is determined based on mouse
% pointer location for most functions, 3) A toggle button for activation
% and deactivation in the figure toolbar, 4) ability to start vcw before
% objects are plotted, 5) "proper" closure of the vcw widget, in fcw the q
% button did not exit the keyDown functions such as panning etc. Now the
% quit action deactivates the widget, 6) Uppon activation of the vcw widget
% the plotting and default view manipulation tools and buttons are disabled
% (to avoid interference with vcw), 7) Added "linked" mode by using ALT
% button to alter views for all axes in figure uppon keypress, 8) Altered
% keypress functions and behaviour with SHIFT, also added i to display help
% information for the vcw function.
% 2015/04/20 Added to GIBBON toolbox
% 2015/04/22 Added JavaFrame handling of ALT related mnemonics
% 2015/04/28 Fixed behaviour for repated vcw; commands (only generate a
% single vcw button even if vcw is called multiple times).
% 2015/04/28 Fixed behaviour for figures without axes. I.e. vcw will only
% start if an axis is present.
% 2016/01/13 Added that clipping is turned off
% 2020/05/11 Have vcw trigger start of vcw when present in figure
% 2020/05/11 Fixed bug in colorbar handling for subplots
% 2020/05/11 Added entries to UserData struct for figure, started adding to
% vcw field. 
% 2020/05/11 Fixed "double start" issue when pressing 'v' 
% 2020/05/11 Have vcw hide the axis interactive toolbar for latest MATLAB
% versions
% 2020/05/11 Removed need for global variable
% 2020/05/11 Fixed bug in scrolling
% 2023/05/29 Added linked mouse based view manipulation using space key,
% removed ALT based linking. 
%
% TO DO: 1) Improved handling of colorbars. Currently requires colorbar
% locations to be set to 'manual' for vcw. However this causes the figure
% to rescale/adjust after deactivation/activation of vcw. It would be best
% if the colorbar locations settings could remain constant. 2) Proper
% restoring of all figure properties. Currently defaults are set manually
% which could remove user defined figure features.
% 2) Turn back the clipping type after exciting vcw
%------------------------------------------------------------------------

%% Parse input arguments
switch nargin
    case 0
        hf = gcf;
        buttonOpt = [];
    case 1
        hf=varargin{1};
        buttonOpt = [];
    case 2
        hf=varargin{1};
        buttonOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end

%Check if hf is a figure handle
if ~isa(hf,'matlab.ui.Figure')
    hf = gcf; %Use current figure (opens new if none exist)
end

%Check view profile option (button mappings)
if ~isa(buttonOpt,'cell')
    if isempty(buttonOpt)
        % Get current view profile (uses default if none is set)
        vcw_profile=getViewProfile;        
    elseif isa(buttonOpt,'char')
        vcw_profile=buttonOpt; %Assume it defines a view profile
    end

    %Set button mapping according to profile
    switch vcw_profile
        case {'CAD','default'}
            buttonOpt = {'pan','rot','zoom','zoom'};
        case 'febio'
            buttonOpt = {'rot','pan','zoom','zoom'};
        case 'touchpad'
            buttonOpt = {'rot','zoom','pan','zoom'};
    end
end

%% Initialise button and button/keypress wait
hb = findall(hf,'Tag','FigureToolBar'); % hb = findall(hf,'Type','uitoolbar');

%Check for presence of a vcw button
hp = findobj(hb,'Tag','tBar');
if isempty(hp) %If vcw button is not present create one and wait for key/button press
    
    % Build icon
    s=[ NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
        NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
        NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;...
        NaN,1  ,1  ,1  ,NaN,1  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN;...
        NaN,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN;...
        1,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN;...
        1,1  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,1  ,NaN,NaN,1  ,1  ,NaN;...
        1,1  ,NaN,NaN,NaN,1  ,1  ,NaN,NaN,1  ,1  ,NaN,NaN,NaN,1  ,1  ;...
        1,1  ,NaN,NaN,NaN,NaN,1  ,NaN,NaN,1  ,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,1  ,1  ,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN;...
        NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,NaN,1  ,1  ,1  ,NaN;...
        NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN;...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN;...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN];
    
    S=zeros(16,16,3);
    S(:,:,1)=0.7.*s;
    S(:,:,2)=0.25.*s;
    S(:,:,3)=0.05.*s;
    
    % Create a uipushtool in the toolbar
    hp=uitoggletool(hb,'TooltipString','Activate View Control Widget (or enter v)','CData',S,'Tag','tBar','Separator','on');
    set(hp,'OnCallback',{@start_vcw,{hf,buttonOpt,hp}});
    set(hp,'OffCallback',{@quit_vcw,{hf,buttonOpt,hp}});
    
    %% Add entries/fields to figure UserData structure
    hf.UserData.vcw.pos=[];
    hf.UserData.vcw.buttonHandle=hp; %The vcw button handle
    hf.UserData.vcw.colorbarHandles=[];
    hf.UserData.vcw.colorbarLocSet={};
    hf.UserData.vcw.linkedOn=-1;

    %% Wait for start using key-press
    
    set(hf,'KeyPressFcn', {@keyPress_wait,buttonOpt,hp,hf},'BusyAction','cancel');
    
else 
    % activate
    drawnow; 
    start_vcw([],[],{hf,buttonOpt,hp});
end

%% Outputs
switch nargout
    case 1
        varargout{1}=hp;
end

end

%%

function keyPress_wait(src,eventData,buttonOpt,hp,hf)

cax = overobj2('axes');

if isempty(cax)
    %     cax=gca; %this gets current axis or if none exists creates one
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end

if isempty(cax)
    return
end

%Key actions
switch eventData.Key
    case {'v'} %Start vcw
        hp.State='On'; % start_vcw(hf,buttonOpt,hp);
end

end
%%
function start_vcw(hObject,callbackdata,inputCell)

hf=inputCell{1};
buttonOpt=inputCell{2};
hp=inputCell{3};

% Store current settings
hf.UserData.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.BusyAction=hf.BusyAction;

% checkAxisLimits(hf);

hp.State='On';
set(hp,'TooltipString','Dectivate View Control Widget (or enter v)');

% Clear any visualization modes we might be in
try
    pan(hf, 'off');
    zoom(hf, 'off');
    rotate3d(hf, 'off');
catch
end

%Quick fix for colorbars
H=findobj(hf,'Type','colorbar'); %Handle set for all colorbars in figure;
if ~isempty(H)
    hf.UserData.vcw.colorbarHandles=H;
    hf.UserData.vcw.colorbarLocSet={H.Location};
    colorbarLocSet(hf,'manual'); %Set colorbar locations to manual
end

% Disable Plottools Buttons and Exploration Buttons
initialState.toolbar = findobj(allchild(hf),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn'),...
        uigettool(initialState.toolbar,'Exploration.Rotate'), ...
        uigettool(initialState.toolbar,'Exploration.Pan'),...
        uigettool(initialState.toolbar,'Exploration.ZoomIn'),...
        uigettool(initialState.toolbar,'Exploration.ZoomOut'),...
        ];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

% For each set of axes
cax = get(hf,'CurrentAxes');
if isempty(cax)
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end

h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
if ~isempty(h)
    for hNow = h
        % Set everything to manual
        set(hNow, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', 'CameraPositionMode', 'manual');
                
        %Set DataAspectRatioMode to manual to keep input aspect ratio
        a=hNow.DataAspectRatio;
        hNow.DataAspectRatio=a;
        hNow.DataAspectRatioMode='manual';
        hNow.DataAspectRatio=a;
        
        % Store the camera viewpoint
        axes(hNow); % axis vis3d;
        caxUserDataStruct.defaultView=camview(hNow);
        set(hNow, 'UserData',caxUserDataStruct);
        
        %Turn clipping off
        set(hNow,'Clipping','off');
   
        %Attemp to turn of interactive toolbar for all axes
        try
            hNow.Toolbar.Visible = 'off';
        catch
        end
        
    end
    axes(cax);
else
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end

% Initialize the callbacks
set(hf, 'WindowButtonDownFcn', {@mousedown, {str2func(['vcw_' buttonOpt{1}]), str2func(['vcw_' buttonOpt{2}]), str2func(['vcw_' buttonOpt{3}])},hf}, ...
    'WindowButtonUpFcn', {@mouseup,hf}, ...
    'KeyPressFcn', {@keypress,buttonOpt,hp,hf}, ...
    'WindowScrollWheelFcn', {@scroll, {str2func(['vcw_' buttonOpt{4}]),hf}}, ...
    'BusyAction', 'cancel');
end

%%

function keypress(~, eventData,buttonOpt,hp,hf)

cax = overobj2('axes');
if isempty(cax)
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end
if isempty(cax)
    return;
end

% checkAxisLimits(hf);

step = 1;
if ismember('shift', eventData.Modifier)
    step = -step; %Make negative while shift is down
end

if ismember('control', eventData.Modifier)
    step = step * 4; %Increase speed
end

% Check for linked view option
linkedOn=hf.UserData.vcw.linkedOn;

% Key input options
switch eventData.Key
    case 'space'
        hf.UserData.vcw.linkedOn=hf.UserData.vcw.linkedOn*-1; %Invert sign
    case 'leftarrow' % Pan left
        vcw_pan([], [step 0], cax, hf);
    case 'rightarrow' % Pan right
        vcw_pan([], [-step 0], cax, hf);
    case 'downarrow' % Pan down
        vcw_pan([], [0 step], cax, hf);
    case 'uparrow' % Pan up
        vcw_pan([], [0 -step], cax, hf);
    case 'x' % Rotate around x
        vcw_rot([], [0 step], cax, hf);
    case 'y' % Rotate around y
        vcw_rot([], [step 0], cax, hf);
    case 'z' % Rotate around z
        vcw_rotz([], [0 step], cax, hf);
    case 'm' % Magnify/zoom (positive or negative)
        vcw_zoom([], [0 -step], cax, hf);
    case 't' % top view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,90);
            end
        else
            view(0,90);
        end
    case 'b' % bottom view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,-90);
            end
        else
            view(0,-90);
        end
    case 'f' % front view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(-90,0);
            end
        else
            view(-90,0);
        end
    case 'h' % back view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(90,0);
            end
        else
            view(90,0);
        end
    case 'l' % left view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,0);
            end
        else
            view(0,0);
        end
    case 'r' % right view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(-180,0);
            end
        else
            view(-180,0);
        end
    case '3' % 3D isometric view 1
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(-45,asind(1/sqrt(3)));
            end
        else
            view(-45,asind(1/sqrt(3)));
        end
    case '4' % 3D isometric view 2
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(135,asind(1/sqrt(3)));
            end
        else
            view(135,asind(1/sqrt(3)));
        end
    case 's' % Store current views as default/initial views (return to it using 'd')
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                caxUserDataStruct=get(h,'UserData');
                caxUserDataStruct.defaultView=camview(h);
                set(h,'UserData',caxUserDataStruct);
            end
        else
            caxUserDataStruct=get(cax,'UserData');
            caxUserDataStruct.defaultView=camview(cax);
            set(cax,'UserData',caxUserDataStruct);
        end
    case 'd' % Reset all the axes to default
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                caxUserDataStruct=get(h,'UserData');
                camview(h,caxUserDataStruct.defaultView);
            end
        else
            caxUserDataStruct=get(cax,'UserData');
            camview(cax,caxUserDataStruct.defaultView);
        end
    case 'i' % Show help
        msgText={'Key input options:',...
            '------------------------------------------------------------------------',...
            'Arrow keys = panning',...
            'x = Rotate around viewer x-axis',...
            'y = Rotate around viewer y-axis',...
            'z = Rotate around viewer z-axis',...
            'm = Zoom/magnify',...
            'Hold down SHIFT to change direction of change for key input based rotation/pan/zoom',...
            'Hold down CTRL to use x4 speed of change  for key input based rotation/pan/zoom',...
            'space = switches on/off linked view mode',...
            'f,h,t,b,l,r = Set front, hind, top, bottom, left, or right view respectively',...
            '3 = 3D isometric view 1',...
            '4 = 3D isometric view 2',...
            's = Store current view states as default (return to default using d)',...
            'd = Restore/reset to default view',...
            'v = activate/deactivate vcw mode',...
            'i = Display help information',...
            'The key inputs and mouse scroll work in the axis defined by mouse pointer location (overobj)',...
            'The mouse inputs (other than scroll) work in current axis (e.g. gca) irrespective of point location',...
            'v = toggle vcw widget on/off',...
            '------------------------------------------------------------------------',...
            };
        helpButton = questdlg(msgText,'Help for vcw','OK','OK');
    case 'v' % Quit vcw mode
        hf.UserData.vcw.linkedOn=-1; %Turn linking off again
        quit_vcw([],[],{hf,buttonOpt,hp});
    otherwise
        % quit_vcw(hf,buttonOpt,hp);
end

end

%%
function mousedown(src, eventData, funcs, hf)

% Get the button pressed
% cax = overobj2('axes');

cax = get(hf, 'CurrentAxes');
if isempty(cax)
    return;
end

% Check for linked view option
linkedOn=hf.UserData.vcw.linkedOn;

% checkAxisLimits(hf);

switch get(hf, 'SelectionType')
    case 'extend' % Middle button
        method = funcs{2};
    case 'alt' % Right hand button
        method = funcs{3};
    case 'open' % Double click
        caxUserDataStructCax=cax.UserData;        
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'    
                caxUserDataStruct=h.UserData; 
                camview(h,caxUserDataStruct.defaultView);
            end
        else
            camview(cax,caxUserDataStructCax.defaultView);
        end

        return;
    otherwise
        method = funcs{1};
end

% Set the cursor
switch func2str(method)
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

% Record where the pointer is
hf.UserData.vcw.pos = get(0, 'PointerLocation');

% Set the cursor and callback
set(hf, 'Pointer', 'custom', 'pointershapecdata', shape, 'WindowButtonMotionFcn', {method, cax,hf});

end

%%
function mouseup(src, eventData,hf)
% Clear the cursor and callback
set(hf, 'WindowButtonMotionFcn', '', 'Pointer', 'arrow');
end

%%
function scroll(src, eventData, inputCell)
func=inputCell{1};
hf=inputCell{2};

% Get the axes handle
cax = overobj2('axes');
if isempty(cax)
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end
if isempty(cax)
    return;
end

% Call the scroll function
func([], [0 -10*eventData.VerticalScrollCount], cax, hf);
end

%%

function d = check_vals(s,d,hf)
% Check the inputs to the manipulation methods are valid

if ~isempty(s)
    % Return the mouse pointers displacement
    new_pt = get(0, 'PointerLocation');
    d = hf.UserData.vcw.pos - new_pt;
    hf.UserData.vcw.pos = new_pt;
end
end

%%

% Figure manipulation functions
function vcw_rot(s, d, cax, hf)
linkedOn=hf.UserData.vcw.linkedOn;
d = check_vals(s,d,hf);
try
    % Rotate XYt
    if linkedOn>0
        for h = findobj(hf, 'Type', 'axes', '-depth', 1)'            
            camorbit(h, d(1), d(2), 'camera', [0 0 1]);
        end
    else
        camorbit(cax, d(1), d(2), 'camera', [0 0 1]);
    end
    axes(cax)
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end

end

function vcw_rotz(s, d, cax, hf)
linkedOn=hf.UserData.vcw.linkedOn;
d = check_vals(s,d,hf);
try
    % Rotate Z
    if linkedOn>0
        for h = findobj(hf, 'Type', 'axes', '-depth', 1)'            
            camroll(h, d(2));
        end
    else
        camroll(cax, d(2));
    end
    axes(cax)
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_zoom(s, d, cax, hf)
linkedOn=hf.UserData.vcw.linkedOn;
d = check_vals(s,d,hf);
% Zoom
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2));
try
    if linkedOn>0
        for h = findobj(hf, 'Type', 'axes', '-depth', 1)'            
            camzoom(h, d);
        end
        axes(cax)
    else
        camzoom(cax, d);
    end
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_zoomz(s, d, cax, hf)
linkedOn=hf.UserData.vcw.linkedOn;
d = check_vals(s,d,hf);
% Zoom by moving towards the camera
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2)) - 1;
try
    if linkedOn>0
        for h = findobj(hf, 'Type', 'axes', '-depth', 1)'            
            camdolly(h, 0, 0, d, 'fixtarget', 'camera');
        end
        axes(cax)
    else
        camdolly(cax, 0, 0, d, 'fixtarget', 'camera');
    end
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_pan(s, d, cax, hf)
linkedOn=hf.UserData.vcw.linkedOn;
d = check_vals(s,d,hf);
try
    % Pan
    if linkedOn>0
        for h = findobj(hf, 'Type', 'axes', '-depth', 1)'            
            camdolly(h, d(1), d(2), 0, 'movetarget', 'pixels');
        end
        axes(cax)
    else
        camdolly(cax, d(1), d(2), 0, 'movetarget', 'pixels');
    end
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function colorbarLocSet(hf,locOpt)
figUserDataStruct=get(hf,'UserData'); %Get user data struct
hColorbarSet=figUserDataStruct.vcw.colorbarHandles; %Color bar handles
for q=1:1:numel(hColorbarSet) %Loop over colorbar handles and set location property
    if isa(locOpt,'cell')
        set(hColorbarSet(q),'Location',locOpt{q});
    else
        set(hColorbarSet(q),'Location',locOpt);
    end
end
end

%%

function quit_vcw(hObject,callbackdata,inputCell)

hf=inputCell{1};
buttonOpt=inputCell{2};
hp=inputCell{3};

% Restore colorbar state
if isfield(hf.UserData.vcw,'colorbarLocSet')    
    colorbarLocSet(hf,hf.UserData.vcw.colorbarLocSet);
end

% Enable Plottools Buttons and Exploration Buttons
try
    initialState.toolbar = findobj(allchild(hf),'flat','Type','uitoolbar');
    if ~isempty(initialState.toolbar)
        initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
            uigettool(initialState.toolbar,'Plottools.PlottoolsOn'),...
            uigettool(initialState.toolbar,'Exploration.Rotate'), ...
            uigettool(initialState.toolbar,'Exploration.Pan'),...
            uigettool(initialState.toolbar,'Exploration.ZoomIn'),...
            uigettool(initialState.toolbar,'Exploration.ZoomOut'),...
            ];
        initialState.ptState = get (initialState.ptButtons,'Enable');
        set (initialState.ptButtons,'Enable','on');
    end
catch
end

set(hp,'State','Off');
set(hp,'TooltipString','Activate View Control Widget (or enter v)');

% Restore figure settings except for key press (if empty) to allow for
% reactivation with v key
hf.WindowButtonDownFcn=hf.UserData.WindowButtonDownFcn;
hf.WindowButtonUpFcn=hf.UserData.WindowButtonUpFcn;
if isempty(hf.UserData.KeyPressFcn)
    hf.KeyPressFcn={@keyPress_wait,buttonOpt,hp,hf};%[];%hf.UserData.KeyPressFcn;
else
    hf.KeyPressFcn=hf.UserData.KeyPressFcn;
end
hf.WindowScrollWheelFcn=hf.UserData.WindowScrollWheelFcn;
hf.BusyAction=hf.UserData.BusyAction;
hf.MenuBar='figure';

%Attemp to turn back on the interactive toolbar for all axes
h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
if ~isempty(h)
    for hNow = h        
        try
            hNow.Toolbar.Visible = 'on';
        catch
        end
    end
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
