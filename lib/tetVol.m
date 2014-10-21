function VE=tetVol(E,V)

% function VE=tetVol(E,V)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the tetrahedral elements specified by the
% element matrix E and the vertices V.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

X=V(:,1); Y=V(:,2); Z=V(:,3);
XE=X(E); YE=Y(E); ZE=Z(E);
if size(E,1)==1 %Transpose in this special case
   XE=XE'; YE=YE'; ZE=ZE';
end
A=[XE(:,1) YE(:,1) ZE(:,1)];
B=[XE(:,2) YE(:,2) ZE(:,2)];
C=[XE(:,3) YE(:,3) ZE(:,3)];
D=[XE(:,4) YE(:,4) ZE(:,4)];

VE=abs(dot((A-D),cross((B-D),(C-D),2),2))./6;
