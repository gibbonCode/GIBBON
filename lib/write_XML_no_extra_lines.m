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
%xmlwrite(save_name,XDOC); %MATLAB's xmlwrite extremely slow for v2017a-2017b
xmlwrite_xerces(save_name,XDOC); %Custom XML write function for now

%Import back into cell array
[T]=txtfile2cell(save_name);

%Save to txt file while skipping "empty lines"
cell2txtfile(save_name,T,1);
 
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
