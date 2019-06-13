function [F,V,regionInd]=multiRegionTriMesh2D(regionSpec,pointSpacing,resampleCurveOpt,plotOn)

% function [F,V,regionInd]=multiRegionTriMesh2D(regionSpec,pointSpacing,plotOn)
% ------------------------------------------------------------------------
% This function creates a 2D triangulation for each of the regions
% specified in the variable regionSpec. The mesh aims to obtain a point
% spacing as defined by the input pointSpacing. The triangulation is based
% on a 2D constrained Delaunary triangulation. The constraints are formed
% by the boundary curves inside the regionSpec cell. Large areas, with
% respect to the pointSpacing, will contain a near homogeneous and
% aproximately equilateral triangulations. Other regions (e.g. at
% boundaries and thin/complex shapes) will contain other triangulation
% types. 
% The function output contains the triangular faces in F, the vertices in V
% and the per triangle region indices in regionInd. By setting plotOn to 0
% or 1 plotting can be switched on or off. 
%
% More on the specification of regions:
% Regions are defined as cell entries in the input variable regionSpec for
% instnace region 1 is found in regionSpec{1}. Each region entry is itself
% also a cell array containing all the boundary curves, e.g. for a two
% curve region 1 we would have something like regionSpec{1}={V1,V2} where
% V1 and V2 are the boundary curves. Multiple curves may be given here. The
% first curve should form the outer boundary of the entire region, the
% curves that follow should define holes inside this boundary and the
% space inside them is therefore not meshed. The boundary vertices for
% regions that share boundaries are merged and will share these boundary
% vertices. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/14/08
%------------------------------------------------------------------------

%% PLOT SETTINGS
if plotOn==1         
    fontSize=20;              
    hf1=cFigure;
    title('Smoothened triangulated mesh','FontSize',fontSize);
    xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
    hold on;  
end

%% MESHING REGIONS

%The total vertex, face and color (=region number) matrices
V=[]; F=[]; regionInd=[];
for qRegion=1:1:numel(regionSpec)
    
    %Define region cell
    regionCell=regionSpec{qRegion};
    
    %Meshing region
    [Fs,Vs]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);
        
    %Joining regions    
    F=[F;Fs+size(V,1)]; %Add new faces and fix vertex indices
    V=[V;Vs]; %Add points
    regionInd=[regionInd; qRegion*ones(size(Fs,1),1)]; %Create region index for faces   
end

%% PLOTTING
if plotOn==1  
    figure(hf1);
    gpatch(F,V,regionInd);
    colormap(gjet(numel(regionSpec))); icolorbar;
    axisGeom(gca,fontSize);    
    drawnow;
end

%% REMOVING DOUBLE POINTS
%Removing double points (region curve points may appear multiple times)
[F,V]=mergeVertices(F,V);
 
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
