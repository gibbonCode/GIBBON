function T = GB_transposefX4(Tin)
%T is a  4X4 transformation matrix of the form  [r r r x
%                                                r r r y
%                                                r r r z
%                                                0 0 0 1]
R = Tin(1:3,1:3);
Ttemp = Tin(1:3,4);

R = R';

Tsub = -R * Ttemp;

T = [R(1,1:3) Tsub(1); R(2,1:3) Tsub(2); R(3, 1:3) Tsub(3); 0 0 0 1];

