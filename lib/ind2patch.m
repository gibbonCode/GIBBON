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
%
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
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
