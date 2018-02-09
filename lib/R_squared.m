function [R_sq]=R_squared(yi,fi)  

% function [R_squared]=R_squared()  
% ------------------------------------------------------------------------
% This function calculates R^2 using the data in vector yi and the modelled
% data in the vector fi. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 05/10/2008
% ------------------------------------------------------------------------

y_bar=nanmean(yi);
SS_err=nansum((yi-fi).^2);
SS_tot=nansum((yi-y_bar).^2);
R_sq=1-(SS_err./SS_tot);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
