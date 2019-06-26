function [R_sq]=R_squared(y,yf)  

% function [R_squared]=R_squared(y,yf)  
% ------------------------------------------------------------------------
% This function calculates the coefficient of determination (R-Squared) for
% the (e.g. measurement) data |y| and the data |yf| (e.g. a model fit). 
% NaN entries in either inputs are ignored. 
%
%
% Change log: 
% 2008/10/05
% 2019/06/26 Updated description
% 2019/06/26 Updated handling of NaN entries in the data
% ------------------------------------------------------------------------

%%

L=~isnan(y) & ~isnan(yf); %Logic for non-NaN entries

y_bar=mean(y(L)); %Mean of y
SS_err=sum((y(L)-yf(L)).^2); %Sum of differences Of y with yf
SS_tot=sum((y(L)-y_bar).^2); %Sum of differences of y with mean
R_sq=1-(SS_err./SS_tot); %R-squared value
 
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
