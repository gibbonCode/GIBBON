function write_XML_no_extra_lines(save_name,XDOC)

% XML_string = xmlwrite(XDOC); %XML as string
% 
% XML_string = regexprep(XML_string,'\n[ \t\n]*\n','\n'); %removes extra tabs, spaces and extra lines
% 
% %Write to file
% fid = fopen(save_name,'w');
% fprintf(fid,'%s\n',XML_string);
% fclose(fid); 

%Write to text file
xmlwrite(save_name,XDOC);

%Import back into cell array
[T]=txtfile2cell(save_name);

%Save to txt file while skipping "empty lines"
cell2txtfile(save_name,T,1);


