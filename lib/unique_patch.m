function [F_uni,V_uni,C_uni,IND_V,IND_F,F_count]=unique_patch(F,V,C,numDigitKeep)

numFacesIni=size(F,1);

%Removing unused vertices
[F,V]=removeNotIndexed(F,V);

%Removing double vertices
try
    [~,IND_V,IND_IND]=unique(round(V,numDigitKeep,'significant'),'rows');
catch
    [~,IND_V,IND_IND]=unique(sround(V,numDigitKeep),'rows');
end
V_uni=V(IND_V,:);
F=IND_IND(F); %Fix indices in F

%Removing double FACES
[F_uni,IND_F,IND_F_2]=uniqueIntegerRow(F);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Fixing face colors, shared faces now obtain mean colour
if ~isempty(C)
    sharedColourMatrixSparse=sparse(IND_F_2,1:numFacesIni,C,numFacesUni,numFacesIni,numFacesIni);    
    C_uni=full(sum(sharedColourMatrixSparse,2))./F_count;
else
    C_uni=[];
end
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
