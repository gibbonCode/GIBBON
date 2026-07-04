function [L]=isGlobalSurfDirOutward(F,V)

% function [L]=isGlobalSurfDirOutward(F,V)
%-------------------------------------------------------------------------
% This function returns a boolean denoting wether surface normals are
% pointing outward (1) (which would result in a positive volume being computed)
% or inward (2) (which would result in a negative volume being computed). 
% The method assumes that the normal directions are coherent across the
% surface. 
%
% Kevin Mattheus Moerman
%
% Change log:
% 2023/08/31 KMM: Updated to use patchVolume
%-------------------------------------------------------------------------
%%

volFV = patchVolume(F,V,0);
L = volFV>0; %Check for positive volume

%% OLD (bad-ish) approach
% Not the best implementation at present. Vertices are offset allong the
% local normal direction by 1/10th of the smalles edge length. Then the
% volume before and after this operation is calculated. If the volume
% decreased the normals face the wrong way and the face orientation is thus
% flipped for smoothening. Contraction/inflation allong normal directions
% in this way does not always yield valid surfaces and hence volume
% computation may be inappropriate.

% %Compute edge lengths
% [edgeLengths]=patchEdgeLengths(F,V);
% growSize=mean(edgeLengths)/10;
% 
% %Get vertex normals
% [~,~,N]=patchNormal(F,V);
% 
% %Compute volumes before and after "contraction/inflation" and flip faces if required.
% [volFV1]=patchVolume(F,V); %Initial volume
% 
% [volFV2]=patchVolume(F,V+(growSize/10.*N)); %"contracted/inflated" volume
% 
% L=volFV2>volFV1; %if the volume increased the global direction is outward

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
