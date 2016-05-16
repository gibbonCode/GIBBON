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
