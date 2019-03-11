function D=import_FBG(loadName)
%-------------------------------------------------------------------------
%
%function D=import_FBG(loadName)
%
%A simple function to import FBG data from text file (specified by
%loadName) and output the data fields (columns) in a structure array D.
%
%
%Kevin Mattheus Moerman
%kevinmoerman@hotmail.com
%2014/01/07 - Created function
%-------------------------------------------------------------------------

%% PARSE TEXT FILE

fid=fopen(loadName);
[cell_out] = textscan(fid,'%s %f %f %f %f %f %f\n', 'delimiter', ',');
fclose(fid);

%% FORMULATE OUTPUT

D.test_date=cell_out{1};
D.cycle_no=cell_out{2};
D.FBG_temp=cell_out{3};
D.FBG_strain=cell_out{4};
D.ind_speed=cell_out{5};
D.ind_depth=cell_out{6};
D.TTL_logic=cell_out{7};

 
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
