function [E,V,CE,Fb,Cb]=hexMeshCubeSphere(cPar)


%% Creating a hollow hexahedral mesh sphere 
% Creating a solid hexahedral mesh sphere

%Control settings
cParShere.sphereRadius=cPar.InnerSphereRadius;
cParShere.coreRadius=cPar.CoreSphereRadius;
cParShere.numElementsMantel=cPar.numElementsSphereCore; 
cParShere.numElementsCore=cPar.numElementsCube; 
cParShere.makeHollow=0;
cParShere.cParSmooth.n=0;

%Creating sphere
[meshStruct]=hexMeshSphere(cParShere);

%Access ouput
E1=meshStruct.E; %The elements 
V1=meshStruct.V; %The vertices
Fb1=meshStruct.Fb;
Cb1=meshStruct.faceBoundaryMarkerBox;

%% Create outer sphere surface and hexahedral mesh


[Fb1c,Vb1c]=patchCleanUnused(Fb1,V1);

[T,P,R] = cart2sph(Vb1c(:,1),Vb1c(:,2),Vb1c(:,3));
[V_crust(:,1),V_crust(:,2),V_crust(:,3)] = sph2cart(T,P,cPar.OuterSphereRadius.*ones(size(R)));

[E2,V2,F2b1,F2b2]=loftLinQuad2hex(Fb1c,Vb1c,V_crust,cPar.numElementsSphereMantel);

%% Create outer box surface mesh

boxDim=cPar.boxWidth*ones(1,3); %Width in each direction
boxEl=cPar.numElementsCube*ones(1,3); %Number of elements per direction 
[meshStruct]=hexMeshBox(boxDim,boxEl);

V_box=meshStruct.V;
F_box=meshStruct.Fb;

%% Make sure the face lists correspond properly

%Remove unused nodes in these lists
[F_box,V_box]=patchCleanUnused(F_box,V_box);
[F_ball,V_ball]=patchCleanUnused(Fb1c,V_crust);

indBall=F_ball(:); 
indBox=F_box(:); 
[~,ind1,~]=unique(indBall); 
indBoxUni=indBox(ind1);
[~,indSort]=sort(indBoxUni);

V_ball=V_ball(indSort,:);
F_ball=indBoxUni(F_ball); 

%% Create hexahedral mesh between outer cube and sphere

[E3,V3,F3b1,F3b2]=loftLinQuad2hex(F_ball,V_ball,V_box,cPar.numElementsCubeSphere);

%% Join and merge node sets

[E,V,CE]=joinElementSets({E1,E2,E3},{V1,V2,V3},{1*ones(size(E1,1),1),2*ones(size(E2,1),1),3*ones(size(E3,1),1)});
[Fb,~,Cb]=joinElementSets({Fb1,F2b2,F3b2},{V1,V2,V3},{7*ones(size(F2b1,1),1),8*ones(size(Fb1,1),1),Cb1});

[E,V,~,ind2]=mergeVertices(E,V,5);
Fb=ind2(Fb);

[FTs,~]=element2patch(E,[],'hex8');

%% Smoothing nodes while holding on to boundaries

cParSmooth.Method='LAP';
cParSmooth.LambdaSmooth=0.5;
cParSmooth.n=cPar.nSmooth;
cParSmooth.RigidConstraints=unique(Fb(:));

if cParSmooth.n>0
    [V]=tesSmooth(FTs,V,[],cParSmooth);
end
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
