function [varargout]=pdfView(varargin)

% function [hFig]=pdfView(fileName,viewerOpt)
%------------------------------------------------------------------------
% View PDF files using a brower embedded in a figure window or within
% MATLAB brower. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/02 : Created based on xmlView
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        fileName=varargin{1};
        viewerOpt=1; 
    case 2
        fileName=varargin{1};
        viewerOpt=varargin{2};
end

%% Prepare html file

pathName=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','temp');
fileName_html=fullfile(pathName,'temp.html');

T={['<embed src="',fileName,'" width="100%" height="100%" alt="pdf" type=''application/pdf''>']};

cell2txtfile(fileName_html,T,0,0);

%%

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
        browser.setCurrentLocation(fileName_html);
        
        %Embed browser
        posPanel = getpixelposition(hFig,true);
        [~,browserContainer] = javacomponent(browser,[1,1,max(posPanel(3)-1,1),max(posPanel(4)-1,1)],hFig);
        set(browserContainer,'Units','normalized');
        drawnow;
    case 2 %Browser viewer
        [~,hFig]=web(fileName_html);%,'-new');
end

if nargout==1
    varargout{1}=hFig;
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
