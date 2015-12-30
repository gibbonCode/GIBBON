function [T]=txtfile2cell(fileName)
% function [T]=txtfile2cell(fileName)
% ------------------------------------------------------------------------
% This function read the text file specified by fileName (path and name)
% whereby each line is read into a seperate entry in the output cell array
% T. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/102/10
%------------------------------------------------------------------------

fid=fopen(fileName);
T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
T=T{1,1};
fclose(fid);