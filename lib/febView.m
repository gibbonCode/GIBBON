function varargout=febView(fileName)


%% View using browser for now

[~,hFig]=web(fileName);%,'-new');

varargout{1}=hFig;

%% Below are test bits to open in figure window in the future

% % T=txtfile2cell(fileName);
% 
% url = fileName;
% 
% hFig = figuremax;
% options.showToolbar = 1; 
% options.showAddressBox = 1; 
% % browser = com.mathworks.mde.webbrowser.WebBrowser.createBrowser(options.showToolbar,options.showAddressBox); 
% browser = com.mathworks.mlwidgets.html.HTMLBrowserPanel;
% % browser.setHtmlText([T{:}]); 
% % pause(0.1); 
% % 
% browser.setCurrentLocation(url);
% posPanel = getpixelposition(hFig,true);
% [~,browserContainer] = javacomponent(browser,[1,1,max(posPanel(3)-1,1),max(posPanel(4)-1,1)],hFig);
% set(browserContainer,'Units','normalized');
% drawnow; 
% varargout{1}=hFig;
% 
% % % Create a figure with a scrollable JEditorPane
% % hfig = figure();
% % je = javax.swing.JEditorPane('text/html', [T{:}]);
% % jp = javax.swing.JScrollPane(je);
% % [hcomponent, hcontainer] = javacomponent(jp, [], hfig);
% % set(hcontainer, 'units', 'normalized', 'position', [0,0,1,1]);
%  
% % % Turn anti-aliasing on (R2006a, Java 5.0)
% % java.lang.System.setProperty('awt.useSystemAAFontSettings', 'on');
% % je.setFont(java.awt.Font('Arial', java.awt.Font.PLAIN, 13));
% % je.putClientProperty(javax.swing.JEditorPane.HONOR_DISPLAY_PROPERTIES, true);
% %  
% % % This only works on Java 1.5 (Matlab R14SP2 to R2007a):
% % je.putClientProperty(com.sun.java.swing.SwingUtilities2.AA_TEXT_PROPERTY_KEY, true);
