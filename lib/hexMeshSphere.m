function [meshStruct]=hexMeshSphere(cPar)

%%

%% Parse input
sphereRadius=cPar.sphereRadius;
coreRadius=cPar.coreRadius;
numElementsMantel=cPar.numElementsMantel;
numElementsCore=cPar.numElementsCore;
makeHollow=cPar.makeHollow;

%% CREATING A HEXAHEDRAL MESH CUBE

%Create box 1
sphereDim=1/sqrt(3)*2*coreRadius*ones(1,3); 
sphereEl=numElementsCore*ones(1,3); %Number of elements
[box2]=hexMeshBox(sphereDim,sphereEl);
E_core=box2.E;
V_core=box2.V;
Fb=box2.Fb;

indBoundary=unique(Fb(:));
V_core_boundary=V_core(indBoundary,:);

%% MAPPING OUTER SURFACE TO A SPHERE

[azimuth,elevation,r] = cart2sph(V_core_boundary(:,1),V_core_boundary(:,2),V_core_boundary(:,3));
[V_core_boundary(:,1),V_core_boundary(:,2),V_core_boundary(:,3)] = sph2cart(azimuth,elevation,coreRadius.*ones(size(r)));
V_core(indBoundary,:)=V_core_boundary;

%% Adding mantel

mantelThickness=sphereRadius-coreRadius;
[Fq,Vq,~]=patchCleanUnused(Fb,V_core);
[E_mantel,V_mantel,F_mantel_inner,F_mantel_outer]=quadThick(Fq,Vq,1,mantelThickness,numElementsMantel);

%Fix outer radii
indBoundary=unique(F_mantel_outer(:));
V_mantel_boundary=V_mantel(indBoundary,:);

[azimuth,elevation,r] = cart2sph(V_mantel_boundary(:,1),V_mantel_boundary(:,2),V_mantel_boundary(:,3));
[V_mantel_boundary(:,1),V_mantel_boundary(:,2),V_mantel_boundary(:,3)] = sph2cart(azimuth,elevation,sphereRadius.*ones(size(r)));
V_mantel(indBoundary,:)=V_mantel_boundary;

%%
if makeHollow==1
    ET=E_mantel; 
    VT=V_mantel; 
    FTb=[F_mantel_outer; F_mantel_inner];
    faceBoundaryMarker=[ones(size(F_mantel_outer,1),1); 2*ones(size(F_mantel_inner,1),1)];
elseif makeHollow==0    
    % Merging node sets
    VT=[V_core;V_mantel];
    ET=[E_core;E_mantel+size(V_core,1)];
    [~,ind1,ind2]=unique(pround(VT,5),'rows');
    VT=VT(ind1,:);
    ET=ind2(ET);
    FTb=ind2(F_mantel_outer+size(V_core,1));
    faceBoundaryMarker=ones(size(F_mantel_outer,1),1);
end

%Get faces
[FT,~]=element2patch(ET,[],'hex8');

%% Smoothing

if ~isfield(cPar,'cParSmooth')
    cPar.cParSmooth.Method='LAP';
    cPar.cParSmooth.LambdaSmooth=0.5;
    cPar.cParSmooth.n=5;
end
cPar.cParSmooth.RigidConstraints=unique(FTb(:));

[F,~,~]=uniqueIntegerRow(FT);
[VT]=tesSmooth(F,VT,[],cPar.cParSmooth);

%% Collect output

meshStruct.E=ET; 
meshStruct.V=VT; 
meshStruct.F=FT;
meshStruct.Fb=FTb;
meshStruct.faceBoundaryMarker=faceBoundaryMarker;

end
