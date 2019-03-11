function [subIndMask,patchDataOut]=sparseSphericalMask(rSphere,nSubTri,optPatchOut)


%Creating icosahedron
[Vp,Fp]=platonic_solid(4,rSphere);

%Subtriangulating
if nSubTri>0
    for q=1:1:nSubTri
        [Fp,Vp]=subtri(Fp,Vp,1);
        [thetaP,phiP,radP] = cart2sph(Vp(:,1),Vp(:,2),Vp(:,3));
        [Vp(:,1),Vp(:,2),Vp(:,3)] = sph2cart(thetaP,phiP,rSphere*ones(size(radP)));
    end
end

subIndMask=unique(round(Vp),'rows');

if optPatchOut==1
    size_M=3*(rSphere)*ones(1,3);
    subIndmaskMiddle=round(size_M*0.5);
    subIndMask=subIndMask+ones(size(Vp,1),1)*(subIndmaskMiddle);
    indMask = unique(sub2ind(size_M, subIndMask(:,2), subIndMask(:,1), subIndMask(:,3)));
    [Fs,Vs,~]=ind2patch(indMask,zeros(size_M),'vu');
    Vs=Vs-ones(size(Vs,1),1)*(subIndmaskMiddle);
    
    patchDataOut.FV_cellTesselation={Fp,Vp};
    patchDataOut.FV_cellMask={Fs,Vs};
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
