function [F,V,F_uv,ij_M]=subdevTexture(F,V,F_uv,ij_M,ns,subdevMethod)

% function [F,V,F_uv,ij_M]=subdevTexture(F,V,F_uv,ij_M,ns,subdevMethod)
% ------------------------------------------------------------------------
% 
% 2023/06/02 KMM: Created for OBJ texture mapping to patch model data
% ------------------------------------------------------------------------

%%

for q=1:1:ns
    switch subdevMethod
        case 1 %Linear
            [Fs1,Vs1]=subtri(F{1},V,1);
            [Fs1_uv,ij_Ms1]=subtri(F_uv{1},ij_M,1);
            [Fs2,Vs2]=subQuad(F{2},V,1);
            [Fs2_uv,ij_Ms2]=subQuad(F_uv{2},ij_M,1);
        case 2 %Loop and Catmull-Clark
            [Fs1,Vs1]=subTriLoop(F{1},V,1,1);
            [Fs1_uv,ij_Ms1]=subTriLoop(F_uv{1},ij_M,1,1);
            [Fs2,Vs2]=subQuadCatmullClark(F{2},V,1,1);
            [Fs2_uv,ij_Ms2]=subQuadCatmullClark(F_uv{2},ij_M,1,1);
    end

    [Fs1,Vs1]=patchCleanUnused(Fs1,Vs1);
    [Fs1_uv,ij_Ms1]=patchCleanUnused(Fs1_uv,ij_Ms1);
    [Fs2,Vs2]=patchCleanUnused(Fs2,Vs2);
    [Fs1_uv,ij_Ms1]=patchCleanUnused(Fs1_uv,ij_Ms1);

    F{1}=Fs1;
    F_uv{1}=Fs1_uv;
    F{2}=Fs2+size(Vs1,1);
    F_uv{2}=Fs2_uv+size(ij_Ms1,1);
    V=[Vs1;Vs2];
    ij_M=[ij_Ms1; ij_Ms2;];
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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
