function [F,V,C]=ind2patch(IND,M,ptype)

% function [F,V,C]=ind2patch(IND,M,ptype)
% ------------------------------------------------------------------------
%
% This function generates patch data (faces 'F', vertices 'V' and color
% data 'C') for 3D images. The patches are only generated for the voxels
% specified by the linear indices in 'IND'. The variable 'ptype' indicates
% the type of patch:
%
% 'v'               Voxel patch data with unshared vertices and faces
%                  such that each voxel has 8 unshared vertices and 6
%                  unshared faces (may be faster than 'vu' and 'vb'
%                  options which require UNIQUE costly computations).
% 'vu'             Voxel patch data such that where possible voxels share
%                  vertices and faces (making patch data computation slower
%                  but plotting more memory efficient).
%                  FaceColor data is averaged for shared faces.
% 'vb'             Voxel patch data faces are exported for unique unshared
%                  faces e.g. only boundary for enclosed volume (plotted
%                  data is visually equivalent to 'v' and 'vu' options when
%                  FaceAlpha is 1)
% 'si', 'sj', 'sk'    Mid-voxel slice patch data for i, j and k direction
%                  respectively
% 'siu', 'sju', 'sku' Same as 'si', 'sj', 'sk' but with double points
%                  removed.
% 'h'              Creates a hexahedral element description instead (e.g
%                  nx8) element data.
% 'hu'             Same as 'h' but with shared unique nodes.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
% 2016/06/10 Added input parsing for IND, i.e. handling empty index set.
% IND will now default to all voxels if input is empty
%------------------------------------------------------------------------

%%

switch ptype
    case 'hu'
        ptype2='h';
    case 'siu'
        ptype2='si';
    case 'sju'
        ptype2='sj';
    case 'sku'
        ptype2='sk';
    otherwise
        ptype2=ptype;
end
[F,V,C]=im2patch(M,IND,ptype2);

%Unshare nodes
switch ptype
    case {'v','h','si','sj','sk'}
        [F,V]=patchDetach(F,V,1);        
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
