function [Ft3]=tri6_subtri3(Ft6,Vt6)

% function [Ft3]=tri6_subtri3(Ft6,Vt6)
% ------------------------------------------------------------------------
% The input 6-node triangles are either oriented such that the first node
% is a corner node or that the first node is a mid-edge node. The two cases
% are highlighted below. 
% 
% Case I
%                     1
%                    / \
%                   2   6
%                  /     \
%                 3___4___5
% Case II
%                     6
%                    / \
%                   5   1
%                  /     \
%                 4___3___2
%
% The output will follow one of these rules and maintains face normals
% directions:
% 
% Case I
%                     1
%                    / \
%                   2___6
%                  / \ / \
%                 3___4___5
%
% Case II
%                     6
%                    / \
%                   5___1
%                  / \ / \
%                 4___3___2
% 
% ------------------------------------------------------------------------
%%
%Check face order based on first
a=Vt6(Ft6(1,2),:)-Vt6(Ft6(1,1),:); %2 1 edge
b=Vt6(Ft6(1,6),:)-Vt6(Ft6(1,1),:); %6 1 edge
c=Vt6(Ft6(1,3),:)-Vt6(Ft6(1,2),:); %3 2 edge

am=sqrt(sum(a.^2));
bm=sqrt(sum(b.^2));
cm=sqrt(sum(c.^2));

phi=acos(dot(a,b)./(am*bm));
phi=mod(phi,pi); %Angle between 2 1 and 6 1 edge 
theta=acos(dot(a,c)./(am*cm));
theta=mod(theta,pi); %Angle between 2 1 and 3 2 edge 

if phi>theta
    Ft3=[Ft6(:,[1 2 6]); Ft6(:,[2 3 4]); Ft6(:,[4 5 6]); Ft6(:,[2 4 6])]; %Clockwise
else 
    Ft3=[Ft6(:,[1 2 3]); Ft6(:,[3 4 5]); Ft6(:,[5 6 1]); Ft6(:,[1 3 5])]; %Anti-clockwise
end


