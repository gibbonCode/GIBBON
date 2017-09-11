function h=axis_inc(s)

%This function widens the axis limits by h_inc=w*(s-1), where w is the
%length of the diagonal of the axis window and s the scaling factor

%%

%Get axis limits
h=axis;

%Calculate increments
wx=h(2)-h(1); 
h_inc_x=wx*(s-1); %Axis limit increment

wy=h(4)-h(3); 
h_inc_y=wy*(s-1); %Axis limit increment

if numel(h)==6
    wz=h(6)-h(5);
    h_inc_z=wz*(s-1); %Axis limit increments
end

%Adjusting extrema
h(1)=h(1)-h_inc_x;
h(2)=h(2)+h_inc_x;
h(3)=h(3)-h_inc_y;
h(4)=h(4)+h_inc_y;
if numel(h)==6
    h(5)=h(5)-h_inc_z;
    h(6)=h(6)+h_inc_z;
end

%Setting axis limits
axis(h);

end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
