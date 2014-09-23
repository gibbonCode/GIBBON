function [logicEuler]=isConvexPolyhedron(F,V)

[eulerVal]=eulerChar(F,V);
 
 logicEuler=eulerVal==2;