%% efw
% Below is a demonstration of the features of the |efw| function

%% Syntax
% |efw(hf);|

%% Description
% The |efw| function, the export figure widget, adds a push button to a
% figure toolbar to link with the |export_fig| function. Press the button
% to start exporting a figure. Users can specify file names, formats,
% resolution and also additional export_fig options. 
%% 
% The Export Figure Widget requires the external function |export_fig|
% created by Oliver Woodford and Yair Altman. It can be obtained from the
% Mathworks Central File Exchange: 
%%
% <http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig>

%% Examples

%% Using the preloaded Export Figure Widget (efw) in |cFigure|
% A |cFigure| window has the Export Figure Widget loaded by default. To use
% |efw| click on the icon in the toolbar. 
% 
% <<efw_icon_menubar.png>>

%%
% Pressing the |efw| button will open a basic |inputdlg| allowing users to
% specify all the usual export_fig options. Hints are given in brackets
% behind the input labels. The default entries are altered according to the
% previous usage within the current figure. 

%%
% <<efw_inputdlg.png>>

%% Use with |figure| 
% The default MATLAB |figure| does not contain the Export Figure Widget. To
% load it here the user must enter |efw;| either before or after axes are
% created. e.g.: 
%%
% |figure; surf(peaks(25)); axis equal; axis tight; efw;|
%%
% which is equivalent to: 
%%
% |figure; efw; surf(peaks(25)); axis equal; axis tight;| 

%% Tip
% You can create your own figure function that simply contains |figure;
% efw;| to create a standard MATLAB figure containing the Export Figure
% Widget.  

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>