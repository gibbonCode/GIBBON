function [varargout]=importFaceFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
faceID=double(A{1});
F=nan(max(faceID),3);
F(faceID,:)=double([A{2} A{3} A{4}]);
faceBoundaryID=double(A{5});

if all(isnan(faceBoundaryID(:)))
    faceBoundaryID=-ones(size(F,1),1);
end

varargout{1}=faceID;
varargout{2}=F;
varargout{3}=faceBoundaryID;

