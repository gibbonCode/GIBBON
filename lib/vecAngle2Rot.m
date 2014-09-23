function [R]=vecAngle2Rot(theta,w)

W=crossProdMat(w);
R=eye(3,3)+W*sin(theta)+W^2*(1-cos(theta));