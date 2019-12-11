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
    % The patch normals are derived by taking the cross product of two edge
    % vectors. Some edge vectors may be (nearly or exactly) co-linear. Hence to
    % select the most suitable edge pair for each face, the edge set with the
    % angle closest to pi/2 is chosen.
    
    %Get indices for most suitable origin and vector pair
    A=abs(patchEdgeAngles(F,V)); %Get edge angles
    A_diff=abs(A-pi/2);%Deviation from 90 degrees
    [~,j]=min(A_diff,[],2); %Get column indices for most oppropriate origin
    JJ=[j(:) j(:)-1 j(:)+1]; %Create previous and next point column indices
    JJ(JJ<1)=size(F,2); %Wrap indices that are too small to end
    JJ(JJ>size(F,2))=1; %Wrap indices that are too large to start
    i=(1:1:size(F,1))'; %Create row index set
    I=i(:,ones(1,size(JJ,2))); %Expand to match the size of the column index set
    indF=sub2ind(size(F),I,JJ); %Convert subscript indices to linear indices in face array
    indSet=F(indF); %Get vertex indices out of face array
    indOrigin=indSet(:,1); %Vertex indices for most suitable origins
    ind1=indSet(:,2); %Vertex indices for first edge point
    ind2=indSet(:,3); %Vertex indices for second edge point
    
    %Compute cross product to get outward normal
    N=cross(V(ind2,:)-V(indOrigin,:),V(ind1,:)-V(indOrigin,:),2);
    N=vecnormalize(N); %Normalize
    
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
