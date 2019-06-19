function gpublish(docName)

%%
if isempty(docName)    
    error('File name is empty');
end

%%

[pathName,docNameClean,~]=fileparts(docName);

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
helpPath=fullfile(toolboxPath,'docs');
htmlPath=fullfile(helpPath,'html');

if isempty(pathName)
    docNameFull=fullfile(helpPath,docName);
else
    docNameFull=docName;
    docName=docNameClean;
end

publish(docNameFull,...
    'catchError',true(1,1),...
    'figureSnapMethod','getframe',...    
    'imageFormat','jpg',...
    'maxHeight',2000,...
    'maxWidth',2000);%

htmlName=fullfile(htmlPath,[docName,'.html']);

[T]=txtfile2cell(htmlName);

imgLineCheck=strfind(T,'<img');

indImgLine=find(~cellfun(@isempty,imgLineCheck));

strAdd=' width="100%" height="auto"';

for q=1:1:numel(indImgLine)   
   lineIndNow=indImgLine(q);
   txtLineNow=T{lineIndNow};
   imgTagLoc=imgLineCheck{lineIndNow};
   indOffset=0;
   for qs=1:1:numel(imgTagLoc)
       indStart=imgTagLoc(qs)+indOffset;
       imgLineCheck_end=regexp(txtLineNow(indStart:end),'>');
       indEnd=imgLineCheck_end(1);
       imgPart=txtLineNow(indStart:indStart+indEnd);       
       if gcontains(imgPart,docName)           
           imgPartNew=[imgPart(1:4),strAdd,imgPart(5:end)];     
           txtLineNow=[txtLineNow(1:indStart-1),imgPartNew,txtLineNow(indStart+indEnd+1:end)];
           indOffset=indOffset+numel(strAdd);
       end       
   end
   T{lineIndNow}=txtLineNow;
      
end
cell2txtfile(htmlName,T,0,0);

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
