function [varargout]=xmlView(fileName)

% function [hFig]=xmlView(fileName)
%------------------------------------------------------------------------
% View XML files using a brower embedded in a figure window
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/05/13 %Updated, renamed from febView to xmlView
%------------------------------------------------------------------------

%%

viewerOpt=1; %Default viewing in figure window
switch viewerOpt
    case 1 %View in figure window
        %Open figure
        figStruct.Name=fileName;
        hFig = cFigure(figStruct);
        
        %Remove tool and menu bars
        ht = findobj(allchild(hFig),'flat','Type','uitoolbar');
        delete(ht);
        ht = findobj(allchild(hFig),'flat','Type','uimenu');
        delete(ht);
        
        %Create browser
        browser = com.mathworks.mlwidgets.html.HTMLBrowserPanel;        
        browser.setCurrentLocation(fileName);
        
        %Embed browser
        posPanel = getpixelposition(hFig,true);
        [~,browserContainer] = javacomponent(browser,[1,1,max(posPanel(3)-1,1),max(posPanel(4)-1,1)],hFig);
        set(browserContainer,'Units','normalized');
        drawnow;
    case 2 %Broswer viewer
        [~,hFig]=web(fileName);%,'-new');
end

if nargout>1
    varargout{1}=hFig;
end
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
