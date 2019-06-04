function [varargout]=textView(varargin)

%%

defaultOptionStruct.BackgroundColor=[44 48 55]/255;
defaultOptionStruct.ForegroundColor=[168 177 190]/255;
defaultOptionStruct.HorizontalAlignment='Left';
defaultOptionStruct.FontName='Monospaced';
defaultOptionStruct.FontWeight='Bold';
defaultOptionStruct.FontSize=12;

switch nargin
    case 1
        T=varargin{1};
        optionStruct=[];
    case 2
        T=varargin{1};
        optionStruct=varargin{2};
end

fileName='';
if ~iscell(T)
    if exist(T,'file')==2
        fileName=T;
        T=txtfile2cell(T);
    end
end

%Fix option structure, complete and remove empty values
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%%

t=(1:1:numel(T));
n=numel(sprintf('%d',max(t)));
t=sprintf(['%-',num2str(n),'.0d    \n'],t); 
t=t(1:end-1);
t=strsplit(t,'\n')';
T=strcat(t,T);

%Open figure
figStruct.Name=fileName;
figStruct.Color=optionStruct.BackgroundColor;
figStruct.MenuBar='none';
figStruct.vcw=0;
figStruct.efw=0;

hf = cFigure(figStruct);

hPan = uipanel(hf, 'Title',fileName, ...
    'Units','pixels','Position',[10 10 hf.Position(3)-20 hf.Position(4)-20],...
    'BorderType','none','BackgroundColor',defaultOptionStruct.BackgroundColor,...
    'ForegroundColor',defaultOptionStruct.ForegroundColor);
hPan.Units='Normalized';
hEdit = uicontrol(hPan, 'Style','edit', 'FontSize',12, ...
    'Min',0, 'Max',2, 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 0 1 1], ...    
    'String',T);

optionSet=fieldnames(optionStruct); 
for q=1:1:numel(optionSet)
   hEdit.(optionSet{q})=optionStruct.(optionSet{q}); 
end

% enable horizontal scrolling
try
jEdit = findjobj(hEdit);
jEditbox = jEdit.getViewport().getComponent(0);
jEditbox.setWrapping(false);                % turn off word-wrapping
jEditbox.setEditable(false);                % non-editable
set(jEdit,'HorizontalScrollBarPolicy',30);  % HORIZONTAL_SCROLLBAR_AS_NEEDED

% maintain horizontal scrollbar policy which reverts back on component resize
hjEdit = handle(jEdit,'CallbackProperties');
set(hjEdit, 'ComponentResizedCallback',...
    'set(gcbo,''HorizontalScrollBarPolicy'',30)')
catch
end
drawnow; 

%%
if nargout>0
    varargout{1}=hf;
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
