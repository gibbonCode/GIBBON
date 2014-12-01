function exportFEB_XML(save_name,XDOC)


% %%
% %Write to text file using xmlwrite (adds extra lines for unknown reasons)
% xmlwrite(save_name,XDOC);
% 
% %Import back into cell array
% [T]=txtfile2cell(save_name);
% 
% %Now save to txt file while skipping "empty lines"
% cell2txtfile(save_name,T,1);

write_XML_no_extra_lines(save_name,XDOC);



