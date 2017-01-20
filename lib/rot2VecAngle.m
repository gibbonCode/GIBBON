function [theta,w]=rot2VecAngle(R)

theta=acos(0.5*(trace(R)-1));
w=(1./(2*sin(theta)))*([R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]);
w=vecnormalize(w);