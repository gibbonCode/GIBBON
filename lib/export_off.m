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
cell2txtfile(fileName,T,0,0);   
 
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
