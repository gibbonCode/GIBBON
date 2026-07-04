function [varargout]=patchEdgeAngles(F,V)

% function [A,Ad]=patchEdgeAngles(F,V)
% -----------------------------------------------------------------------
% Computes the edge angles (A) for the patch data specified by the faces
% (F) and vertices (V) arrays.
%
%
% See also:
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/08/18
%------------------------------------------------------------------------

%%

numNodes=size(F,2);

A=zeros(size(F));
for q=1:1:numNodes
    
    q1=q;
    q0=q-1;
    if q0<1
        q0=numNodes;
    end
    q2=q+1;
    if q2>numNodes
        q2=1;
    end
    
    P0=V(F(:,q0),:)-V(F(:,q1),:);
    P0=vecnormalize(P0);
    P2=V(F(:,q2),:)-V(F(:,q1),:);
    P2=vecnormalize(P2);
    A(:,q)=acos(dot(P0,P2,2));
    
end

varargout{1}=A;

if nargout>1
    varargout{2}=180*(A./pi);
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
