function [R T] = fX4_to_RT(Tm)

R = Tm(1:3,1:3);
T = Tm(1:3,4)';