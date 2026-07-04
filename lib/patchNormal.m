function [varargout]=patchNormal(F,V)

% [N,Vn,Nv]=patchNormal(F,V)
% ------------------------------------------------------------------------
% Normals are derived based on cross product of triangle edge vectors. Each
% triangle is constructed using the two points of its face and the mean of
% the face.
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/06/02 Updated general lay-out and commenting
% 2015/09/22 Fixed to allow for 2D patch data
% 2016/11/15 Added handling of non-triangular faces
% 2019/04/22 Change to use cross product of edges closest to orthogonal
% 2019/04/23 Added handling of cell input
% 2021/09/13 Updated to use more accurate (for non-planar patches) and efficient patchEdgeCrossProduct method
%------------------------------------------------------------------------

%%
% Deal with cell input by calling function recursively
if isa(F,'cell')
    N=repmat({zeros(size(F,1),3)},size(F));
    for q=1:1:numel(F)
        N{q}=patchNormal(F{q},V);
    end
    Vn=patchCentre(F,V); %Face center cell
    if nargout==3
        Nv=faceToVertexMeasure(F,V,N); %Normal vectors at vertices
    end
else
    
    %%
    % Deal with 2D patch data
    if size(V,2)==2
        V(:,3)=0;
    end

    %% 
    %Compute patch face normal direction    
    
    %Compute edge vector cross products
    C=patchEdgeCrossProduct(F,V); %0.5 sum of edge cross products
    
    %Get normal vectors through normalization
    N=vecnormalize(C); %Normalization
    
    %%
    
    %Calculate face centres if requested
    if nargout>1
        %Get mean face coordinates for normal vectors
        Vn=patchCentre(F,V);
    end
    
    %Calculate vertex normals if requested
    if nargout==3
        [Nv]=faceToVertexMeasure(F,V,N); %Convert face data to vertex data
        Nv=vecnormalize(Nv); %Normalize vectors
    end
        
end

%% Collect output

switch nargout
    case 1 %Only face normals
        varargout{1}=N;
    case 2 %Face normals and face centres
        varargout{1}=N;
        varargout{2}=Vn;
    case 3 %Face normals, face centres and vertex normals
        varargout{1}=N;
        varargout{2}=Vn;
        varargout{3}=Nv;
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
