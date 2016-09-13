%% cell2txtfile
% Below is a demonstration of the features of the |cell2txtfile| function

%% Syntax
% |cell2txtfile(fileName,T,skipOpt);|

%% Description
% The |cell2txtfile| function exports the content in the cell array T to
% the text file fileName. Each entry in the cell array will be a line in
% the txt file. Prior to text file creation the cell is converted to a
% column format. If the input skipOpt=1 cell entries which appear empty
% (after spaces are removed) will be skipped.

%%
clear; close all; clc;

%% Examples

%% Exporting a cell containing text data to a txt file  
% Create example cell containing text entries, empty entries, an entry
% containing just spaces and a non-char entry number. The number entry will
% be converted to a string. 

T={'Hello','','World','  ',uint8(125)};
filePath=mfilename('fullpath');
fileName=fullfile(fileparts(filePath),'data','temp','temp.txt');
skipOpt=0; %Empty entries will be kept 
cell2txtfile(fileName,T,skipOpt);

%%
%Output text reread as cell:
[T_out]=txtfile2cell(fileName);
disp(T_out)

%%
% Changing skipOpt to 1 will skip these empty lines
skipOpt=1; %Empty entries will be skipped
cell2txtfile(fileName,T,skipOpt);

[T_out]=txtfile2cell(fileName);
disp(T_out)

%%
% Remove example file
delete(fileName);

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
