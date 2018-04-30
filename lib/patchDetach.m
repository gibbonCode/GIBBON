function [Fc,Vc]=patchDetach(F,V,shrinkFactor)

%%

if numel(shrinkFactor)==size(V,1) %If specified on the nodes
    shrinkFactor=vertexToFaceMeasure(F,shrinkFactor); % Convert to face metric
end

if numel(shrinkFactor)~=1 && numel(shrinkFactor)~=size(F,1)
    error('The number of elements in shrinkFactor should be equal to 1 or the number of faces or the number of nodes');
end

%%
Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    FX_mean=mean(FX,2);
    FX=((FX-FX_mean).*shrinkFactor)+FX_mean;
    Vc(:,q)=FX(:);
end
    
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
