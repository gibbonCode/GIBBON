function [D]=patchEdgeLengths(F,V)

% function [D]=patchEdgeLengths(F,V)
% -----------------------------------------------------------------------
% Computers the edge lengths (D) for the patch data specified by the faces
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
% kevinmoerman@hotmail.com
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

%Computer the edge lengths
D=sqrt(sum(VD.^2,2));

