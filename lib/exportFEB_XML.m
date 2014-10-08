function exportFEB_XML(save_name,XDOC)

%Write to text file using xmlwrite (adds extra lines for unknown reasons)
xmlwrite(save_name,XDOC);

%Import back into cell array
[T]=txtfile2cell(save_name);

%Add GIBBON comment at start
TT=cell(1,1);
TT(1,1)=T(1,1);
TT(2,1)={'<!-- '};
TT(3,1)={['Created using GIBBON, ',datestr(now)]};
TT(4,1)={'-->'};
TT(end+1:end+size(T,1)-1)=T(2:end);

%Now save to txt file while skipping "empty lines"
cell2txtfile(save_name,TT,1);


