function [Ft,Vt]=hemiSphereCylMesh(S)

% function [Ft,Vt]=hemiSphereCylMesh(S)
% ------------------------------------------------------------------------
% Creates the faces (Ft) and vertices (Vt) for a triangular surface
% describing a hemisphere connected to a cylinder. The input is a structure
% array S containint the following control parameters:
%   S.sphereRadius => The radius of the hemi-spher portion
%   S.nRefineRegions => Number of |subtri| refinements for icosahedron
%   S.cylinderHeight => height of the cylinder part
%   S.cylinderStepSize => Aproximate node spacing for cylinder portion
%
% See also: |hemiSphereRegionMesh|
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/13/08
% ------------------------------------------------------------------------

%% Defining hemi-sphere portion
patchType=S.patchType;
switch patchType
    case 'tri' %Triangular sphere mesh
        hemiStruct.sphereRadius=S.sphereRadius; %Sphere radius
        hemiStruct.nRefineRegions=S.nRefine; %Number of refinement steps for regions
        hemiStruct.nRefineMesh=1; %Number of refinement steps for mesh
        [Fs,Vs,~]=hemiSphereRegionMesh(hemiStruct);
        Fs=fliplr(Fs); %flip face orientation
        Ft=Fs;
    case 'quad' %Create quadrilateral sphere mesh        
        [Vs,Fs]=platonic_solid(2,S.sphereRadius);
        Fs=fliplr(Fs); %flip face orientation
        for q=1:1:S.nRefine
            [Fs,Vs]=subQuad(Fs,Vs,1);
            [T,P,R] = cart2sph(Vs(:,1),Vs(:,2),Vs(:,3));
            [Vs(:,1),Vs(:,2),Vs(:,3)]=sph2cart(T,P,ones(size(R)).*S.sphereRadius);
        end
        
        %Get top hemi-sphere
        LV=Vs(:,3)>(0-eps(0));
        LF=all(LV(Fs),2);        
        Fs=Fs(LF,:); %Crop faces
        
        %Remove unused vertices and fix face indices accordingly
        uniInd=unique(Fs(:));
        indAll=nan(1,size(Vs,1));
        indAll(uniInd)=1:numel(uniInd);
        Vs=Vs(uniInd,:);
        Fs=indAll(Fs);
        Ft=[Fs(:,1) Fs(:,3) Fs(:,2); Fs(:,1) Fs(:,4) Fs(:,3)];
        Ft=fliplr(Ft); %flip face orientation
end
Vs(:,3)=Vs(:,3)-max(Vs(:,3)); %Set max at zero
Vs(:,3)=-Vs(:,3); %Flip around

%% Defining cylindrical portion

%Cylinder parameters
cylinderHeight=S.cylinderHeight; %Cylinder height
cylinderStepSize=S.cylinderStepSize; %Edge length in Z direction for cylinder

% Find hemi-sphere edge
TR =triangulation(Ft,Vs);
edgeDome = freeBoundary(TR);
indEdges=unique(edgeDome(:));
V_edge=Vs(indEdges,:);

%Reorder edge points to create curve
% [V_edge,~]=curvePathOrderFix(V_edge);
[x,y,~]=getColumns(V_edge); %Get coordinates
T = atan2(y,x); %polar angles 
[~,indSort]=sort(T); %Sort angles
V_edge=V_edge(indSort,:);
[x,y,z]=getColumns(V_edge); %Reorder and get coordinates

%Determine Z step sizes
if isempty(cylinderStepSize);
    cylinderStepSize=mean(sqrt(sum(diff(V_edge,1,1).^2,2))); %mean point spacing on edge
end
nStepsCylinder=round(cylinderHeight./cylinderStepSize);

%Create cylinder mesh from edge points
X=x(:,ones(1,nStepsCylinder));
Y=y(:,ones(1,nStepsCylinder));
Z=linspacen(z,z+cylinderHeight,nStepsCylinder);

X=X'; Y=Y'; Z=Z';
%Create patch data
[F,Vc] = surf2patch(X,Y,Z);

%Stitch ends together
I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
F_sub1=sub2ind(size(Z),I,J);
Fq=[F;F_sub1];

%%
switch patchType
    case 'tri' %Triangular sphere mesh
        %Convert to triangles
        Fc=[Fq(:,1) Fq(:,3) Fq(:,2); Fq(:,1) Fq(:,4) Fq(:,3)];
        Fc=fliplr(Fc); %flip face orientation
    case 'quad' %Create quadrilateral sphere mesh
        Fc=Fq;
end

%% MERGING MODEL PORTIONS

%Combining vertices
Vt=[Vs;Vc];

%Combining faces
Fcc=Fc+size(Vs,1); %Fix vertex ID's here
Ft=[Fs;Fcc];

%Removing double vertices and fixing face indices
[~,ind1,ind2]=unique(round(Vt*1e5),'rows');
Vt=Vt(ind1,:);
Ft=ind2(Ft);