function [out1,out2,out3] = qginput(varargin)

% [out1,out2,out3] = qginput(varargin)
% ------------------------------------------------------------------------
% Modified version of |ginput| function. The difference is that this
% function lets you specify the mouse pointer type and by default avoids the
% use of the 'fullcrosshair' mouse pointer since it seems to cause errors
% when combined with patch graphics. Another difference is the exit key.
% For MATLAB's ginput the Enter key will exit ginput. For qginput this has
% been replaced by the Esc key. 
% varargin{1} is assumed similar to the input n for ginput while
% varargin{2} if present will be the mouse pointer type. If not present
% then the mousePointer variable will be set to the default: 'smallHand'.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/10
% 2016/06/27 KMM Changed "break" char to Esc rather than Enter
%------------------------------------------------------------------------

%% 
% JUNK CODE:
% %Could replace ginput by something allong these lines, however there
% seems to be more to it than this i.e. in relation to whether computer is
% a PC or not and what the figure units are. 

% %Allocate output vectors
% x=zeros(numEvents,1);
% y=zeros(numEvents,1);
% k=zeros(numEvents,1);
% for q=1:1:numEvents
%
%     %Alter pointer
%     set(gcf,'pointer','cross');
%     drawnow;
%
%     %Wait for click/key stroke
%    [keyDown]=waitForValidButtonPress;
%    k(q)=keyDown;
%
%    %Get coordinates
%    P = get(gca, 'CurrentPoint');
%
%    x(q)=
%
% end

%% PARSE INPUT

switch nargin
    case 1
        arg1=varargin{1};
%         mousePointerType='crosshair'; %Default
        [mousePointerType]=specialPointerShape('smallHand'); %Default
    case 2
        arg1=varargin{1};
        mousePointerType=varargin{2}; %User specified mouse pointer (empty results in no mouse pointer alteration)
    otherwise
        error(message('MATLAB: False number of inputs', 'qginput'))
end

%% REST IS SIMILAR TO GINPUT
out1 = []; out2 = []; out3 = []; y = [];

if ~matlab.ui.internal.isFigureShowEnabled
    error(message('MATLAB:hg:NoDisplayNoFigureSupport', 'qginput'))
end

fig = gcf;
figure(gcf);

if nargin == 0
    how_many = -1;
    b = [];
else
    how_many = arg1;
    b = [];
    if  ischar(how_many) ...
            || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
            || ~(fix(how_many) == how_many) ...
            || how_many < 0
        error(message('MATLAB:qginput:NeedPositiveInt'))
    end
    if how_many == 0
        % If input argument is equal to zero points,
        % give a warning and return empty for the outputs.
        
        warning (message('MATLAB:qginput:InputArgumentZero'));
    end
end

% Setup the figure to disable interactive modes and activate pointers.
initialState = setupFcn(fig,mousePointerType); %KMM added mousePointerType input

% onCleanup object to restore everything to original state in event of
% completion, closing of figure errors or ctrl+c.
c = onCleanup(@() restoreFcn(initialState));


% We need to pump the event queue on unix
% before calling WAITFORBUTTONPRESS
drawnow
char = 0;

while how_many ~= 0
    % Use no-side effect WAITFORBUTTONPRESS
    waserr = 0;
    try
        keydown = wfbp;
    catch %#ok<CTCH>
        waserr = 1;
    end
    if(waserr == 1)
        if(ishghandle(fig))
            cleanup(c);
            error(message('MATLAB:qginput:Interrupted'));
        else
            cleanup(c);
            error(message('MATLAB:qginput:FigureDeletionPause'));
        end
    end
    % g467403 - qginput failed to discern clicks/keypresses on the figure it was
    % registered to operate on and any other open figures whose handle
    % visibility were set to off
    figchildren = allchild(0);
    if ~isempty(figchildren)
        ptr_fig = figchildren(1);
    else
        error(message('MATLAB:qginput:FigureUnavailable'));
    end
    %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
    %         clicked figure has handlevisibility set to callback
    if(ptr_fig == fig)
        if keydown
            char = get(fig, 'CurrentCharacter');
            button = abs(get(fig, 'CurrentCharacter'));
        else
            button = get(fig, 'SelectionType');
            if strcmp(button,'open')
                button = 1;
            elseif strcmp(button,'normal')
                button = 1;
            elseif strcmp(button,'extend')
                button = 2;
            elseif strcmp(button,'alt')
                button = 3;
            else
                error(message('MATLAB:qginput:InvalidSelection'))
            end
        end
        axes_handle = gca;
        drawnow;
        pt = get(axes_handle, 'CurrentPoint');
        
        how_many = how_many - 1;

        %% KMM Changed "break" char to Esc rather than Enter
        if(char == 27) % & how_many ~= 0)
            % if the ESC key was pressed, char will ==27,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.
            % If this was an early breakout, don't include
            % the <Esc> key info in the return arrays.
            % We will no longer count it if it's the last input.            
            break; 
        end
        %% OLD VERSION
%         if(char == 13) % & how_many ~= 0)
%             % if the return key was pressed, char will == 13,
%             % and that's our signal to break out of here whether
%             % or not we have collected all the requested data
%             % points.
%             % If this was an early breakout, don't include
%             % the <Return> key info in the return arrays.
%             % We will no longer count it if it's the last input.
%             
%             break; 
%         end
%%
        
        out1 = [out1;pt(1,1)]; %#ok<AGROW>
        y = [y;pt(1,2)]; %#ok<AGROW>
        b = [b;button]; %#ok<AGROW>
    end
end

% Cleanup and Restore
cleanup(c);

if nargout > 1
    out2 = y;
    if nargout > 2
        out3 = b;
    end
else
    out1 = [out1 y];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok<NASGU>

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error(message('MATLAB:qginput:Interrupted'));
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function initialState = setupFcn(fig,mousePointerType)

% Store Figure Handle.
initialState.figureHandle = fig;

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

% Setup mouse pointer without warning.
if ~isempty(mousePointerType)
    oldwarnstate = warning('off', 'MATLAB:Figure:Pointer');
    
    %Set user specified mousePointerType
    if ischar(mousePointerType) %normal options
         set(fig,'Pointer',mousePointerType); 
    else %its assumed that the input is a structure specifying the pointer
        
    %mousePointerType.PointerShapeCData;%A 16x16 array defining a custom pointer CDATA
    %mousePointerType.PointerShapeHotSpot;%Row column indices of custom pointer centre    
    set(fig,'Pointer','custom','PointerShapeCData',mousePointerType.PointerShapeCData,'PointerShapeHotSpot',mousePointerType.PointerShapeHotSpot); %User specified mousePointerType
    
    end
    
   
    warning(oldwarnstate);
end

% Adding this to enable automatic updating of currentpoint on the figure
set(fig,'WindowButtonMotionFcn',@(o,e) dummy());

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
end

function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    set(initialState.figureHandle,'WindowButtonMotionFcn','');
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
end

function dummy()
% do nothing, this is there to update the GINPUT WindowButtonMotionFcn.
end

function cleanup(c)
if isvalid(c)
    delete(c);
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
