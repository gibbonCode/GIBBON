function [DI]=TriScatteredInterp_ND(DT,D,XI,InterpMethod)

if ~isa(DT,'DelaunayTri') %if DT is not a delaunay tesselation
    DT=delaunayTriangulation(DT); %assuming DT are coordinates, replace by Delaunay tesselation
end

DI=nan(size(XI,1),size(D,2)); %Allocate DI
for q=1:size(D,2)% loop over dimensions
    switch InterpMethod
        case 'nat_near' %natural in chull, neirest outside chull
            [DI(:,q),~]=TriScatteredInterp_nat_near(DT,D(:,q),XI);
        otherwise %TriScatterInterp can handle other methods
            F = scatteredInterpolant(DT,D(:,q),InterpMethod);
            %Construct interpolator
            DI(:,q)=F(XI); %Interpolate
    end
end

end