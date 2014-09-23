function [r,c]=thetaphi2vec(THETA,PHI)

%--------------------------------------------------------------------------
% function [vx,vy]=thetaphi2vec(THETA,PHI)
% 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/02/2010
%-------------------------------------------------------------------------- 


%%

r = [cos(THETA).*cos(PHI)  -sin(THETA) cos(THETA).*sin(PHI)];
c = [sin(THETA).*cos(PHI)  cos(THETA)  sin(THETA).*sin(PHI)];

end
