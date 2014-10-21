function [varargout]=importNodeFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %f %f %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
nodeID=double(A{1});
V=nan(max(nodeID),3);
V(nodeID,:)=[A{2} A{3} A{4}];

varargout{1}=nodeID;
varargout{2}=V;

