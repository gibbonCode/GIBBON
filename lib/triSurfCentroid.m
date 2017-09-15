function [Vc]=triSurfCentroid(F,V)

%%
A=patch_area(F,V);
Vcc=patchCentre(F,V);
Vc=sum(Vcc.*A,1)./sum(A(:)); 
