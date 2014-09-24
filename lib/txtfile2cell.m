function [T]=txtfile2cell(filename)

fid=fopen(filename);

if ~isempty(strfind(version,'R2014a'))
    T=textscan(fid,'%s','delimiter', '\n','Whitespace','','bufsize',1e6);
else
    T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
end

T=T{1,1};
fclose(fid);