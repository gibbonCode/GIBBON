function [V_MP,Dc]=inddisp2cumdisp(V,Di,chull_crop,interp_method)


%24/10/2013 Changed Delaunay methods due to MATLAB version changes

%%

n_ini=size(V,3);
n_disp=size(Di,3);

%% Creating Delaunay tesselations of initials

L_valid_INI=cell(1,n_ini);
DT_INI=cell(1,n_ini);
for q=1:1:n_ini
    X=V(:,:,q);
    L_valid_INI{q}=~any(isnan(X),2);
    DT_INI{q} = delaunayTriangulation(X(L_valid_INI{q},:));
end

if n_ini==1; %If only one initial is provided we replicate the initials and derived parameters to the size of the displacement array
    L_valid_INI=repmat(L_valid_INI,[1 n_disp]);
    DT_INI=repmat(DT_INI,[1 n_disp]);
    V=repmat(V,[1 1 n_disp]);
end

%% Deriving cumulative displacement through interpolation

%Allocating memory
Dc=NaN(size(Di));
Dc(:,:,1:2)=Di(:,:,1:2); %First dynamic is initial, second is first valid displacement
V_MP=NaN(size(V));
V_MP(:,:,1)=V(:,:,1);%Initial coordinates
V_MP(:,:,2)=V(:,:,1)+Di(:,:,2); %First material point set

for i_dyn=3:n_disp;
     
    XI=V_MP(:,:,i_dyn-1); %Previous material point locations    
    X=V(:,:,i_dyn); %Current initial coordinates 
    DX=Di(:,:,i_dyn); %Current displacement measure from initial
    
    L_valid_1=~any(isnan(DX),2);
    L_valid_2=~any(isnan(XI),2);
    if chull_crop==1
        %Check if these points are within the convex hull of the initial coordinates
%         L_inchull=~isnan(tsearchn(DT_INI{i_dyn}.Points,DT_INI{i_dyn}.ConnectivityList,XI(L_valid_2,:)));
        L_inchull = ~isnan(pointLocation(DT_INI{i_dyn},XI(L_valid_2,:)));
    else
        L_inchull=true(size(XI(L_valid_2,1)));
    end    
    LXI_inchull=false(size(XI)); LXI_inchull(L_valid_2)=L_inchull; %Convert back to logic index in XI
    L_valid_3=L_valid_INI{i_dyn} & L_valid_1; %Logic for valid initial points AND valid displacements
    
    if strcmp(interp_method,'v4'); %biharmonic spline based
%         DX_m=biharm_spline_fit_ND(X(L_valid_3,:),DX(L_valid_3,1),XI(LXI_inchull,:));
%         DY_m=biharm_spline_fit_ND(X(L_valid_3,:),DX(L_valid_3,2),XI(LXI_inchull,:));
%         DZ_m=biharm_spline_fit_ND(X(L_valid_3,:),DX(L_valid_3,3),XI(LXI_inchull,:));
%         DI_inchull=[DX_m DY_m DZ_m];
    elseif strcmp(interp_method,'nat_near');
        %         DT=delaunayTriangulation(X(L_valid_3,:));
        %         [Dx,Lnn]=TriScatteredInterp_nat_near(X(L_valid_3,:),DX(L_valid_3,1),XI(LXI_inchull,:));
        %         [Dy,Lnn]=TriScatteredInterp_nat_near(X(L_valid_3,:),DX(L_valid_3,2),XI(LXI_inchull,:));
        %         [Dz,Lnn]=TriScatteredInterp_nat_near(X(L_valid_3,:),DX(L_valid_3,3),XI(LXI_inchull,:));
        %         DI_inchull=[Dx Dy Dz];
        FX = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,1),'natural','nearest');
        FY = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,2),'natural','nearest');
        FZ = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,3),'natural','nearest');
        DI_inchull=[FX(XI(LXI_inchull,:)) FY(XI(LXI_inchull,:)) FZ(XI(LXI_inchull,:))];
    else%Delaunay based interpolation
        %Construct interpolators for X,Y and Z
%         DT=delaunayTriangulation(X(L_valid_3,:));
        FX = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,1),interp_method,'none');
        FY = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,2),interp_method,'none');
        FZ = scatteredInterpolant(X(L_valid_3,:),DX(L_valid_3,3),interp_method,'none');                
        DI_inchull=[FX(XI(LXI_inchull,:)) FY(XI(LXI_inchull,:)) FZ(XI(LXI_inchull,:))];        
    end
    
    DI=NaN(size(XI)); DI(LXI_inchull,:)=DI_inchull; 
    V_MP(:,:,i_dyn)=V_MP(:,:,i_dyn-1)+DI;
    Dc(:,:,i_dyn)=DI;
    
end

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
