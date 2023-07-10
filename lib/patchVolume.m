function [surfaceVolume]=patchVolume(varargin)
% function [surfaceVolume]=patchVolume(F,V,makeAbsolute)
% ------------------------------------------------------------------------
% Volume derivation based on Gauss divergence theorem
%
%
%
%
% ------------------------------------------------------------------------
%% parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        makeAbsolute=0;
    case 3
        F=varargin{1};
        V=varargin{2};
        makeAbsolute=varargin{3};
end

%%

if isa(F,'cell')
    surfaceVolume=0;
    for q=1:1:numel(F)
        if isa(V,'cell')
            [surfaceVolumeContributions]=getVolumeContributions(F{q},V{q},makeAbsolute);
        else
            [surfaceVolumeContributions]=getVolumeContributions(F{q},V,makeAbsolute);
        end
        %Total volume. Ignore NaNs, which may result from normal computation on zero-area faces
        surfaceVolume = surfaceVolume + sum(surfaceVolumeContributions(~isnan(surfaceVolumeContributions)));
    end
else
    [surfaceVolumeContributions]=getVolumeContributions(F,V,makeAbsolute);
end
%Total volume. Ignore NaNs, which may result from normal computation on zero-area faces
surfaceVolume = sum(surfaceVolumeContributions(~isnan(surfaceVolumeContributions))); 

end

function [surfaceVolumeContributions]=getVolumeContributions(F,V,makeAbsolute)    
    N=patchNormal(F,V);%Face normals
    surfaceAreas=patchArea(F,V); %Face areas
    Z=V(:,3);
    Zm=mean(Z(F),2); %Mean Z-coordinates for faces            
    Nz = N(:,3); %Z component of normal
    surfaceVolumeContributions = surfaceAreas.*Zm.*Nz; %Contributions
    if makeAbsolute
        surfaceVolumeContributions=abs(surfaceVolumeContributions);
    end
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