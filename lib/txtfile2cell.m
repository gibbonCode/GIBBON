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
% 2014/10/10
%------------------------------------------------------------------------

fid=fopen(fileName);

if ~isempty(strfind(version,'R2014a'))
    T=textscan(fid,'%s','delimiter', '\n','Whitespace','','bufsize',1e6);
else
    T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
end

T=T{1,1};
fclose(fid);