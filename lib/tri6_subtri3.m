function [Ft3]=tri6_subtri3(Ft6)

Ft3=[Ft6(:,[1 2 6]); Ft6(:,[2 3 4]); Ft6(:,[4 5 6]); Ft6(:,[2 4 6])];