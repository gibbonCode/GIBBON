function [W]=mask_design(Z)

% function [Im Jm Km Wm]=mask_design(I, J, K, Z)
% ------------------------------------------------------------------------
% This function designs the weights for a mask. It assumes the vector 'Z'
% contains the normalised intensities found in a region of interest. The
% mask weight vector 'W' is calculated such that the value sum(W.*Z) is
% minimum when the mask is centered on the region of interest and maximum
% when it is not.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 13/08/2008
% ------------------------------------------------------------------------


%%

W=Z;
W(Z==1)=-(10./(numel(Z(Z==1))*(Z(Z==1))));
W((Z>0 & Z<1))=(10./(numel(Z(Z>0 & Z<1))*(Z((Z>0 & Z<1)))));
W(Z==0)=999;

%% END
 
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
