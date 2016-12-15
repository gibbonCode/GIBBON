function export_off(fileName,F,V)


%%
%Get edges
[E]=patchEdges(F,1);

%Create top text
T={'OFF ';...
       [num2str(size(V,1)),' ',num2str(size(F,1)),' ',num2str(size(E,1))];...
       ''};
   
%Create vertex text
textForm=repmat('%f ',1,size(V,2));
textForm=[textForm,'\n'];
TV=sprintf(textForm,V');
T{end+1}=TV(1:end-1); %Take off added end of line statement

%Create faces text
F_mat=[size(F,2)*ones(size(F,1),1) F-1]; 
textForm=repmat('%u ',1,size(F,2)+1);
textForm=[textForm,'\n'];
TF=sprintf(textForm,F_mat');
T{end+1}=TF;

%Write text to file
cell2txtfile(fileName,T,0);   