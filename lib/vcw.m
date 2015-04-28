function vcw(varargin)

% function vcw(hf,buttonOpt)
% ------------------------------------------------------------------------
% vcw, View Control Widget, Allows users to manipulate a view in 3D using key and button presses
%
%   vcw(hf,buttonOpt)
%
% Allows the user to rotate, pan and zoom a figure using key presses and
% mouse gestures. Additionally, press q to quit the widget, r to reset the
% axes and escape to close the figure. This function is non-blocking, but
% fixes axes aspect ratios.
%
% IN:
%   hf - Handle of the figure to be manipulated (default: gcf).
%   buttonOpt - 4x1 cell array indicating the function to associate with
%             each mouse button (left to right) and the scroll action.
%             Functions can be any of:
%                'rot' - Rotate about x and y axes of viewer's coordinate
%                        frame
%                'rotz' - Rotate about z axis of viewer's coordinate frame
%                'zoom' - Zoom (change canera view angle)
%                'zoomz' - Move along z axis of viewer's coordinate frame
%                'pan' - Pan
%                '' - Don't use that button
%             Default: {'pan','rot','zoomz','zoomz'};).
%
% This code was inspired by the fcw function by Oliver Woodford (which was
% based on Torsten Vogel's view3d function, which was in turn inspired by
% rotate3d from The MathWorks, Inc.).
%
% The vcw function differs from the fcw function in the following:
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/04/15 %Copied from fcw and renamed to vcw due to planning revision
% 2015/04/15 %Added: 1) handing of colorbars (bug in fcw when view(2) is
% used combined with panning which induced zooming and panning), 2) overobj
% axes selection so that the current axes is determined based on mouse
% pointer location for most functions, 3) A toggle button for activation
% and deactivation in the figure toolbar, 4) ability to start vcw before
% objects are plotted, 5) "proper" closure of the vcw widget, in fcw the q
% button did not exit the keyDown functions such as panning etc. Now the
% quit action deactivates the widget and waits for reactivation, 6) Uppon
% activation of the vcw widget the plotting and default view manipulation
% tools and buttons are disabled (to avoid interference with vcw), 7) Added
% "linked" mode by using ALT button to alter views for all axes in figure
% uppon keypress, 8) Altered keypress functions and behaviour with SHIFT,
% also added i to display help information for the vcw function.
% 2015/04/20 Added to GIBBON toolbox
% 2015/04/22 Added JavaFrame handling of ALT related mnemonics
% 2015/04/28 Fixed behaviour for repated vcw; commands (only generate a
% single vcw button even if vcw is called multiple times). 
% 2015/04/28 Fixed behaviour for figures without axes. I.e. vcw will only
% start if an axis is present. 
%
% TO DO: 1) Improved handling of colorbars. Currently requires colorbar
% locations to be set to 'manual' for vcw. However this causes the figure
% to rescale/adjust after deactivation/activation of vcw. It would be best
% if the colorbar locations settings could remain constant. 2) Proper
% restoring of all figure properties. Currently defaults are set manually
% which could remove user defined figure features.
%------------------------------------------------------------------------

%% Parse input arguments
switch nargin
    case 0
        hf = gcf;
        buttonOpt = {'pan','rot','zoomz','zoomz'};
    case 1
        hf=varargin{1};
        buttonOpt = {'pan','rot','zoomz','zoomz'};
    case 2
        hf=varargin{1};
        buttonOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end

if ~ishandle(hf)
    buttonOpt=hf;
    hf = gcf;
end


%% Initialise button and button/keypress wait
hb = findall(hf,'Type','uitoolbar');

%Check for presence of a vcw button
hp = findobj(hb,'Tag','tBar');

if isempty(hp); %If vcw button is not present create one and wiat for key/button press
    
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
    hp=uitoggletool(hb,'TooltipString','Activate View Control Widget (or enter v)','CData',S,'Tag','tBar');
    set(hp,'OnCallback',{@start_vcw_toggle,{hf,buttonOpt,hp}});
    set(hp,'OffCallback',{@quit_vcw_toggle,{hf,buttonOpt,hp}});
    
    %% Wait for start using key-press
    
    set(hf,'KeyPressFcn', {@keyPress_wait,buttonOpt,hp},'BusyAction','cancel');
    
end
return

function start_vcw_toggle(hObject,callbackdata,x)
start_vcw(x{1},x{2},x{3});
return

function keyPress_wait(src,eventData,buttonOpt,hp)
hf = ancestor(src, 'figure');
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
        start_vcw(hf,buttonOpt,hp);
end
return

function start_vcw(hf,buttonOpt,hp)

set(hp,'State','On');
set(hp,'TooltipString','Dectivate View Control Widget (or enter v)');

% Clear any visualization modes we might be in
pan(hf, 'off');
zoom(hf, 'off');
rotate3d(hf, 'off');

%Quick fix for colorbars
H=findobj(gcf,'Type','colorbar'); %Handle set
figUserDataStruct=get(hf,'UserData');
if isempty(figUserDataStruct)
    figUserDataStruct.colorbarLocSet=get(H,'Location');
    set(hf,'UserData',figUserDataStruct);
end
colorbarLocSet(gcf,'manual');

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
% cax=gca; %this gets current axis or if none exists creates one
if isempty(cax)
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end
h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles

if ~isempty(h)
    for h = findobj(hf, 'Type', 'axes', '-depth', 1)';
        % Set everything to manual
        set(h, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', 'CameraPositionMode', 'manual');
        % Store the camera viewpoint
        axes(h); axis vis3d;
        caxUserDataStruct.defaultView=camview(h);
        set(h, 'UserData',caxUserDataStruct);
    end
    axes(cax);
else
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end

% Initialize the callbacks
set(hf, 'WindowButtonDownFcn', {@mousedown, {str2func(['vcw_' buttonOpt{1}]), str2func(['vcw_' buttonOpt{2}]), str2func(['vcw_' buttonOpt{3}])}}, ...
    'WindowButtonUpFcn', @mouseup, ...
    'KeyPressFcn', {@keypress,buttonOpt,hp}, ...
    'WindowScrollWheelFcn', {@scroll, str2func(['vcw_' buttonOpt{4}])}, ...
    'BusyAction', 'cancel');
return

function keypress(src, eventData,buttonOpt,hp)

hf = ancestor(src, 'figure');
cax = overobj2('axes');
if isempty(cax)
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end
if isempty(cax)
    return;
end

step = 1;
if ismember('shift', eventData.Modifier)
    step = -step; %Make negative while shift is down
end

if ismember('control', eventData.Modifier)
    step = step * 4; %Increase speed
end

mnemOff=1;
if ismember('alt', eventData.Modifier)
    linkedOn=1;
    % Try to turn off menu Mnemonics
    try
        warning off; %Stop jframe warning
        jFrame = get(handle(hf),'JavaFrame');
        jMenuBar=jFrame.fHG2Client.getMenuBar;
        for q=0:1:jMenuBar.getComponentCount-1;
            jComp=jMenuBar.getComponent(q);
            jComp.setMnemonic(' ');
        end
        mnemOff=1;
        warning on;
    catch
        %Remove toolbar when ALT is pressed (QUICK FIX)
        mnemOff=0;
        set(hf,'MenuBar','none');
        t = uitoolbar;
        set(t,'Tag','emptyBar_vcw');
    end
else
    linkedOn=0;
end

% Key input options
switch eventData.Key
    case 'leftarrow' % Pan left
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [step 0], h);
            end
        else
            vcw_pan([], [step 0], cax);
        end
    case 'rightarrow' % Pan right
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [-step 0], h);
            end
        else
            vcw_pan([], [-step 0], cax);
        end
    case 'downarrow' % Pan down
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [0 step], h);
            end
        else
            vcw_pan([], [0 step], cax);
        end
    case 'uparrow' % Pan up
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [0 -step], h);
            end
        else
            vcw_pan([], [0 -step], cax);
        end
    case 'x' % Rotate around x
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rot([], [0 step], h);
            end
        else
            vcw_rot([], [0 step], cax);
        end
    case 'y' % Rotate around y
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rot([], [step 0], h);
            end
        else
            vcw_rot([], [step 0], cax);
        end
    case 'z' % Rotate around z
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rotz([], [0 step], h);
            end
        else
            vcw_rotz([], [0 step], cax);
        end
    case 'm' % Magnify/zoom (positive or negative)
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_zoom([], [0 -step],h);
            end
        else
            vcw_zoom([], [0 -step], cax);
        end
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
            'Hold down ALT to link manipulations for all figure axes for key input based rotation/pan/zoom',...
            'f,h,t,b,l,r = Set front, hind, top, bottom, left, or right view respectively',...
            's = Store current view states as default (return to default using d)',...
            'd = Restore/reset to default view',...
            'v = activate/deactivate vcw mode',...
            'q = Deactivate vcw mode and unable/remove widget from current figure',...
            'i = Display help information',...
            'The key inputs and mouse scroll work in the axis defined by mouse pointer location (overobj)',...
            'The mouse inputs (other than scroll) work in current axis (e.g. gca) irrespective of point location',...
            '------------------------------------------------------------------------',...
            };
        helpButton = questdlg(msgText,'Help for vcw','OK','OK');
    case 'v' % Quit vcw mode
        quit_vcw(hf,buttonOpt,hp);
    case 'q' % Quit vcw mode
        quit_vcw(hf,buttonOpt,hp);
        close_vcw(hf,hp);
    case 'escape'
        close(hf);
end
if mnemOff==0
    set(hf,'MenuBar','figure');
    h1 = findobj(hf,'Tag','emptyBar_vcw');
    if ~isempty(h1)
        delete(h1);
    end
end
return

function mousedown(src, eventData, funcs)

% Get the button pressed
hf = ancestor(src, 'figure');
cax = get(hf, 'CurrentAxes');
if isempty(cax)
    return;
end

switch get(hf, 'SelectionType')
    case 'extend' % Middle button
        method = funcs{2};
    case 'alt' % Right hand button
        method = funcs{3};
    case 'open' % Double click
        caxUserDataStruct=get(cax,'UserData');
        camview(cax,caxUserDataStruct.defaultView);
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
global VCW_POS
VCW_POS = get(0, 'PointerLocation');
% Set the cursor and callback
set(ancestor(src, 'figure'), 'Pointer', 'custom', 'pointershapecdata', shape, 'WindowButtonMotionFcn', {method, cax});

return

function mouseup(src, eventData)
% Clear the cursor and callback
set(ancestor(src, 'figure'), 'WindowButtonMotionFcn', '', 'Pointer', 'arrow');
return

function scroll(src, eventData, func)
% Get the axes handle
hf = ancestor(src, 'figure');
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
func([], [0 -10*eventData.VerticalScrollCount], cax);
return

function d = check_vals(s, d)
% Check the inputs to the manipulation methods are valid
global VCW_POS
if ~isempty(s)
    % Return the mouse pointers displacement
    new_pt = get(0, 'PointerLocation');
    d = VCW_POS - new_pt;
    VCW_POS = new_pt;
end
return

% Figure manipulation functions
function vcw_rot(s, d, cax)
d = check_vals(s, d);
try
    % Rotate XY
    camorbit(cax, d(1), d(2), 'camera', [0 0 1]);
catch
    % Error, so release mouse down
    mouseup(cax)
end
return

function vcw_rotz(s, d, cax)
d = check_vals(s, d);
try
    % Rotate Z
    camroll(cax, d(2));
catch
    % Error, so release mouse down
    mouseup(cax)
end
return

function vcw_zoom(s, d, cax)
d = check_vals(s, d);
% Zoom
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2));
try
    camzoom(cax, d);
catch
    % Error, so release mouse down
    mouseup(cax)
end
return

function vcw_zoomz(s, d, cax)
d = check_vals(s, d);
% Zoom by moving towards the camera
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2)) - 1;
try
    camdolly(cax, 0, 0, d, 'fixtarget', 'camera');
catch
    % Error, so release mouse down
    mouseup(cax)
end
return

function vcw_pan(s, d, cax)
d = check_vals(s, d);
try
    % Pan
    camdolly(cax, d(1), d(2), 0, 'movetarget', 'pixels');
catch
    % Error, so release mouse down
    mouseup(cax)
end
return

function colorbarLocSet(hf,locOpt)
H=findobj(hf,'Type','colorbar'); %Handle set
for q=1:1:numel(H)
    if isa(locOpt,'cell')
        set(H(q),'Location',locOpt{q});
    else
        set(H(q),'Location',locOpt);
    end
end
return

function quit_vcw_toggle(hObject,callbackdata,x)
quit_vcw(x{1},x{2},x{3})
return

function quit_vcw(hf,buttonOpt,hp)

% Restore figure settings except for key press to allow reactivation
set(hf, 'WindowButtonDownFcn',[], ...
    'WindowButtonUpFcn',[], ...
    'KeyPressFcn',{@keyPress_wait,buttonOpt,hp}, ...
    'WindowScrollWheelFcn',[], ...
    'BusyAction','cancel',...
    'MenuBar','figure');

% Restore colorbar state
figUserDataStruct=get(hf,'UserData');
colorbarLocSet(hf,figUserDataStruct.colorbarLocSet);

% Enable Plottools Buttons and Exploration Buttons
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

%         hp=findobj(hf,'Tag','tBar');
set(hp,'State','Off');
set(hp,'TooltipString','Activate View Control Widget (or enter v)');

return

function close_vcw(hf,hp)

% Restore figure settings
set(hf, 'WindowButtonDownFcn',[], ...
    'WindowButtonUpFcn',[], ...
    'KeyPressFcn',[], ...
    'WindowScrollWheelFcn',[], ...
    'BusyAction','queue',...
    'MenuBar','figure');

% Restore colorbar state
figUserDataStruct=get(hf,'UserData');
colorbarLocSet(hf,figUserDataStruct.colorbarLocSet);

% Enable Plottools Buttons and Exploration Buttons
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

%         hp=findobj(hf,'Tag','tBar');
set(hp,'State','Off');
set(hp,'TooltipString','Activate View Control Widget (or enter v)');
delete(hp);
return