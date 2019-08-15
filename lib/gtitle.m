function [varargout]=gtitle(varargin)

% function hText=gtitle(titleString,fontSize,hf)
% -----------------------------------------------------------------------
%
%
% Change log: 
% 2019/08/09 Wrapped try statement around JavaFrame based background color
% removal and suppressed warning occuring in MATLAB R2019b
% -----------------------------------------------------------------------
%% Parse input

switch nargin
    case 1
        titleString=varargin{1};
        fontSize=[];
        hf=[];     
        optionStruct=[];
    case 2
        titleString=varargin{1};
        fontSize=varargin{2};
        hf=[];
        optionStruct=[];
    case 3
        titleString=varargin{1};
        fontSize=varargin{2};
        hf=varargin{3};
        optionStruct=[];
    case 4
        titleString=varargin{1};
        fontSize=varargin{2};
        hf=varargin{3};
        optionStruct=varargin{4};
end

%%
defaultOptionStruct.HorizontalAlignment='Center';
defaultOptionStruct.FontWeight='bold';
defaultOptionStruct.FontName='Helvetica';
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%%

if isempty(fontSize)
    fontSize=11;
end

if isempty(hf)
    hf=gcf;
end

backGroundColor=hf.Color;

if strcmp(backGroundColor,'k')
    textColor='w';
elseif size(backGroundColor,2)==3
    grayLevels=linspace(0,1,100);
    [~,indMax]=max(abs(grayLevels-mean(double(backGroundColor))));    
    textColor=grayLevels(indMax)*ones(1,3);
else
    textColor='k';    
end

hText = uicontrol(hf,'Style','text','String',titleString,'BackgroundColor',backGroundColor,...
    'HorizontalAlignment',optionStruct.HorizontalAlignment,'FontSize',fontSize,...
    'FontWeight',optionStruct.FontWeight,'Units','Points','ForegroundColor',textColor,...
    'FontName',optionStruct.FontName);

% Attempt to turn of background color
try %JavaFrame approach
    warning 'off'
    j_hText = findjobj(hText);
    j_hText.setOpaque(false);
    j_hText.repaint();
    warning 'on'
catch
    % Not possible to make background transparent yet without JavaFrame
end

hText.Units = 'Points';

figResize([],[],{hf,hText});

hFunc=get(hf,'ResizeFcn');

if iscell(hFunc)
%     warning('gtitle replaced the ResizeFcn function. Specify your ResizeFcn in the form @(h,e)figResize(h,e,c) to avoid this behavior');    
    set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hText}));
else
    if isempty(hFunc)
        set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hText}));
    else        
        set(hf,'ResizeFcn',@(a,b)(cellfun(@(x)feval(x,a,b),{hFunc,@(a,b)figResize(a,b,{hf,hText})})));
    end
end

%% Collect output

if nargout>0
    varargout{1}=hText;
end

end

function figResize(~,~,inputCell)
hf=inputCell{1};
hTextInfo=inputCell{2};
unitsNow=hf.Units;
hf.Units='Points';

figPosition=hf.Position;
textBoxWidth=figPosition(3);
textBoxHeight=hTextInfo.Extent(4).*ceil(hTextInfo.Extent(3)./figPosition(3));
textPosition=[0 figPosition(4)-textBoxHeight textBoxWidth textBoxHeight];
hTextInfo.Position = textPosition;

hf.Units=unitsNow;
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
