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
