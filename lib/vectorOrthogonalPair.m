function [a,d]=vectorOrthogonalPair(f)

[max_f,J_max]=max(f,[],2); %Largest entries in vector set
J_max=J_max+1; %Add one to the column index
J_max(J_max>3)=J_max(J_max>3)-2; %Go the otherway if out of bounds
I_max=(1:1:size(f,1))'; %Row index
indMax=sub2ind(size(f),I_max,J_max); %Linear index

a=f; %a initialized f 
a(indMax)=a(indMax)+max_f; %a forced different by adding maximum entry one of its directions
[a]=vecnormalize(a); %a normalised and not equal to f

d=cross(f,a); [d]=vecnormalize(d); %d is orthogonal to f and a
a=cross(d,f); [a]=vecnormalize(a); %a is reset to be orthogonal to both f and d
