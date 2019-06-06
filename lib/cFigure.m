function [varargout]=cFigure(varargin)

% function [h]=cFigure(figStruct)
% ------------------------------------------------------------------------
% Creates a custom figure using the input structure figStruct. The cFigure
% function provides easy control of background color, the color
% definitions, the figure window size (e.g. near maximal), and enables
% figure property cloning. It also allows users to create hidden figures
% which can be made visible for instance using the mfv command.
%
% The content of figStruct may follow all properties of a normal figure
% i.e. such that figStruct=figure. Which could lead to (amonst other
% properties):
%
% figStruct=figure
%
%   Figure (1) with properties:
%
%       Number: 10
%         Name: ''
%        Color: [0.9400 0.9400 0.9400]
%     Position: [680 558 560 420]
%        Units: 'pixels'
%
% figStruct used to be a handle in which case its use in this function
% involves the set command e.g.: set(h,'outerPosition',[a b c d]);
% For newer MATLAB versions however the cFigure function uses a different
% but equivalent syntax i.e.:
% h.outerPosition=[a b c d];
%
% Some additional fields can be added that are not normally part of the
% figure property set: ColorDef and ScreenOffset.
% ColorDef sets the color definition which is either 'white' or 'black'.
% This allows the user to select a dark background and appropriately set
% the colorscheme for it e.g. for a black background:
%         figStruct.ColorDef='black';
%         figStruct.Color='k';
% Where the Color property sets the figure background color while the
% ColorDef property sets the colorscheme used (of axes etc.).
% By default cFigure creates figures that are the full screensize but
% reduced 10% away from the edges. The spacing between the figure window
% and the screen edges is set by the figStruct.ScreenOffset property. The
% units are pixels.
%
% See also: figure, set, get, colordef, mfv, scf
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log:
% 2014/11/25 Created
% 2015/04/15 Added vcw functionality
% 2018/02/02 Fixed bug in relation to groot units (e.g. figure size is
% wrong if units are not pixels). 
% 2018/12/05 Updated figure size handling which now accepts either a
% screen offset or a scaling factor. 
% 2018/12/19 Fixed bug to handle when both offset and scaling are specified
% 
% To do: 
% Check handling of multiple screens for figure size setting
%------------------------------------------------------------------------

%% Parse input and set defaults

%Force groot units to be pixels
graphicalRoot=groot;
grootUnits=graphicalRoot.Units;
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units='pixels';
end
screenSizeGroot = graphicalRoot.ScreenSize(3:4); %Get screen widht and height

%Default settings
defaultFigStruct.Visible='on';
defaultFigStruct.ColorDef='white';
defaultFigStruct.Color='w';
defaultFigStruct.ScreenScale=0.85; %Figure size is based on scaled screensize
% defaultFigStruct.ScreenOffset=round(max(screenSizeGroot)*0.1); %i.e. figures are spaced around 10% of the sreensize from the edges        
defaultFigStruct.Clipping='off';
defaultFigStruct.efw=1;
defaultFigStruct.vcw={'pan','rot','zoomz','zoomz'};

switch nargin
    case 0
       figStruct=[]; %Use default
    case 1        
        figStruct=varargin{1}; %Use custom
end

%Fix option structure, complete and remove empty values
[figStruct]=structComplete(figStruct,defaultFigStruct,1);

%Get export figure option and remove field
efwOpt=figStruct.efw;
figStruct=rmfield(figStruct,'efw'); %Remove field from structure array

%Get view control widget options and remove field
vcwOpt=figStruct.vcw;
figStruct=rmfield(figStruct,'vcw'); %Remove field from structure array
            
%%

isOld=verLessThan('matlab', '8.4.0.150421 (R2014b)');

%% Create a hidden figure
hf = figure('Visible', 'off'); %create an invisible figure

%% Setcolor definition and associated defaults
hf=colordef(hf,figStruct.ColorDef); %Update figure handle
figStruct=rmfield(figStruct,'ColorDef'); %Remove field from structure array

%% Set figure size

%Set figure units
if isOld
    set(hf,'units','pixels');    
else
    hf.Units='pixels';    
end

%Set figure size
if isfield(figStruct,'ScreenOffset')        
    %Compute figure offset from border        
    figSizeEdgeOffset=figStruct.ScreenOffset/2; 
elseif isfield(figStruct,'ScreenScale')
    %Compute screen offset
    figSizeEdgeOffset=(screenSizeGroot-(screenSizeGroot*figStruct.ScreenScale))/2;        
end

%Remove custom fields from structure
if isfield(figStruct,'ScreenOffset')        
    figStruct=rmfield(figStruct,'ScreenOffset'); %Remove field from structure array
end

if isfield(figStruct,'ScreenScale')
    figStruct=rmfield(figStruct,'ScreenScale'); %Remove field from structure array
end

%Get/set units
figUnits=hf.Units; %Get current figure units (users may change defaults)
hf.Units=graphicalRoot.Units; %Force units the same

%Compute figure size in terms of width and height
figSize=screenSizeGroot-figSizeEdgeOffset*2; % width, height offsets

%Position and resize figure
if isOld
    set(hf,'outerPosition',[(screenSizeGroot(1)-figSize(1))/2 (screenSizeGroot(2)-figSize(2))/2 figSize(1) figSize(2)]); % left bottom width height
else
%     hf.Position=[figSizeEdgeOffset(1) figSizeEdgeOffset(2) figSize(1) figSize(2)]; % left bottom width height
    hf.Position=[(screenSizeGroot(1)-figSize(1))/2 (screenSizeGroot(2)-figSize(2))/2 figSize(1) figSize(2)]; % left bottom width height
end

%Set figure units back
hf.Units=figUnits;

%% Parse remaining figure properties

% Note: This is where figure becomes visible if figStruct.Visible='on'

fieldSet = fieldnames(figStruct); % Cell containing all structure field names
for q=1:1:numel(fieldSet)
    fieldNameCurrent=fieldSet{q};
    try
        if isOld
            set(hf,fieldNameCurrent,figStruct.(fieldNameCurrent));
        else
            hf.(fieldNameCurrent)=figStruct.(fieldNameCurrent);
        end
    catch errorMsg
        rethrow(errorMsg); %likely false option
    end
end

%% Check for activation of vcw

if isa(vcwOpt,'cell') %Allow enabling of vcw mode    
    hp=vcw(hf,vcwOpt);
    hf.UserData.cFigure.Handles.vcw=hp;
end

%% Check for activation of efw
if efwOpt
    efw; 
end

%%
if nargout>0
    varargout{1}=hf;
end

%%
% Reset groot units if a change was needed
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units=grootUnits;
end

%%
% Initiate a set of drawnow events to cope with MATLAB bug (observed in
% MATLAB 2019a). Figure position might alter after several drawnow events. 

pos=hf.Position;
for q=1:1:10    
    drawnow;    
    if ~all(pos==hf.Position)        
        break
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
