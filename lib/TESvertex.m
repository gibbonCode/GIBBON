function [TES_dist]=TESvertex(TES,x,y,z)

x_vert=x(TES);
y_vert=y(TES);
z_vert=z(TES);

c_from=1:1:size(TES,2);
c_upto=[c_from(end) c_from(1:end-1)];
TES_dist=zeros(size(TES,1),numel(c_from));
for c=c_from;
    dc=sqrt((x_vert(:,c_from(c))-x_vert(:,c_upto(c))).^2+(y_vert(:,c_from(c))-y_vert(:,c_upto(c))).^2+(z_vert(:,c_from(c))-z_vert(:,c_upto(c))).^2);
    TES_dist(:,c)=dc;
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
