function [R,Xc]=insphere(E,X)


TR = triangulation(E,X); %Get triangulated object
Xc = incenter(TR); %Calculate circum center coordinates
V  = tetVol(E,X); %Get element volumes

%Compute element areas
A1 = patch_area(E(:,[1 2 3]),X);
A2 = patch_area(E(:,[1 2 4]),X);
A3 = patch_area(E(:,[2 3 4]),X);
A4 = patch_area(E(:,[3 1 4]),X);
A  = A1 + A2 + A3 + A4; %Element areas

R  = 3*V./A;
