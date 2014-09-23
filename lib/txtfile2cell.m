function [T]=txtfile2cell(filename)

fid=fopen(filename);
T=textscan(fid,'%s','delimiter', '\n','Whitespace','','bufsize',1e6);
T=T{1,1};
fclose(fid);