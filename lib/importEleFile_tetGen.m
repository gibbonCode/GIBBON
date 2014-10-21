function [varargout]=importEleFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
elementID=double(A{1});
E=nan(max(elementID),4);
E(elementID,:)=double([A{2} A{3} A{4} A{5}]);
elementMaterialID=double(A{6});

if all(isnan(elementMaterialID(:)))
    elementMaterialID=-ones(size(E,1),1);
end

varargout{1}=elementID;
varargout{2}=E;
varargout{3}=elementMaterialID;

