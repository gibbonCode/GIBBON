function [L]=isGlobalSurfDirOutward(F,V)

% function [L]=isGlobalSurfDirOutward(F,V)
%-------------------------------------------------------------------------
% Not the best implementation at present. Vertices are offset allong the
% local normal direction by 1/10th of the smalles edge length. Then the
% volume before and after this operation is calculated. If the volume
% decreased the normals face the wrong way and the face orientation is thus
% flipped for smoothening. Contraction/inflation allong normal directions
% in this way does not always yield valid surfaces and hence volume
% computation may be inappropriate.
%-------------------------------------------------------------------------
%%

%Compute edge lengths
[edgeLengths]=patchEdgeLengths(F,V);
minLength=min(edgeLengths(:));

%Get vertex normals
[~,~,N]=patchNormal(F,V);

%Compute volumes before and after "contraction/inflation" and flip faces if required.
[volFV1]=triSurfVolume(F,V); %Initial volume
[volFV2]=triSurfVolume(F,V+(minLength/10.*N)); %"contracted/inflated" volume

L=volFV2>volFV1; %if the volume increased the global direction is outward

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
