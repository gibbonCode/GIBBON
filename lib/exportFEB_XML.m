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
