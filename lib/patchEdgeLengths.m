function [D]=patchEdgeLengths(F,V)

% function [D]=patchEdgeLengths(F,V)
% -----------------------------------------------------------------------
% Computes the edge lengths (D) for the patch data specified by the faces
% (F) and vertices (V) arrays. If size(F,2)>2 it is assumed that F indeed
% represents faces. If however size(F,2)==2 it is instead assumed that F is
% an array representing edges. As such it skips the computation of the
% edges array. The edges array used is non-unique by default. See the
% |patchEdges| function for more details if the lengths of a unique set of
% edges is desired. 
%
%
% See also: |patchEdges|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/17
%------------------------------------------------------------------------

%%

%Derive edge array
if size(F,2)>2 %The input is assumed to represent faces hence an edge array is derived
    E=patchEdges(F);
else %It is assumed that the input array represents an edges array
    E=F; 
end

%Derive edge vertex arrays
V_E1=V(E(:,1),:);
V_E2=V(E(:,2),:);

%Derive difference vectors
VD=(V_E1-V_E2);

%Compute the edge lengths
D=sqrt(sum(VD.^2,2));

 
%% 
% ********** _license boilerplate_ **********
% 
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
