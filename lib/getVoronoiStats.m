function voronoiStats=getVoronoiStats(Fv_cell,Vv,np,numBins)

%Replace infinite entries by NaN
Vv(isinf(Vv))=NaN;

%Allocate memory for structure
voronoiStats.ellipse=nan(numel(Fv_cell),5);
voronoiStats.ellipticity=nan(numel(Fv_cell),1);
voronoiStats.circularity=nan(numel(Fv_cell),1);
voronoiStats.angle=nan(numel(Fv_cell),1);
voronoiStats.area=nan(numel(Fv_cell),1);
voronoiStats.radii=cell(numel(Fv_cell),1);
voronoiStats.minRad=nan(numel(Fv_cell),1);
voronoiStats.maxRad=nan(numel(Fv_cell),1);
voronoiStats.meanRad=nan(numel(Fv_cell),1);

%Loop through cells to compute statistics
for q=1:1:numel(Fv_cell);
    
    %Get current Voronoi cell
    fv=Fv_cell{q};
    vv=Vv(fv,:);
    
    if ~isnan(vv(:))
        %Compute ellips and ellipticity
        [A] = ellipseFit(vv,2,np);
        
        a=A(3); b=A(4);
        if b<a %Oblate
            E=sqrt((a^2-b^2)./a^2);
        else %Prolate
            E=sqrt((b^2-a^2)./b^2);
        end
        
        %Compute radii
        x=vv(:,1);
        y=vv(:,2);
        xCentre = mean(x);
        yCentre = mean(y);
        x = (x-xCentre);
        y = (y-yCentre);
        R=sqrt(x.^2+y.^2);

        %Compute area
        vv_area=patch_area(fv,Vv);

        %Compute circularity
        C=polyCircularity(vv);
        
        %Compute angle
        ellipseAngle=A(end);
        ellipseAngle(ellipseAngle<0)=(2*pi+ellipseAngle(ellipseAngle<0))-pi;
        
        %Store in stucture
        voronoiStats.ellipse(q,:)=A;
        voronoiStats.ellipticity(q)=E;
        voronoiStats.circularity(q)=C;
        voronoiStats.angle(q)=ellipseAngle;
        voronoiStats.area(q)=vv_area;
        voronoiStats.radii{q}=R;
        voronoiStats.minRad(q)=min(R);
        voronoiStats.maxRad(q)=max(R);
        voronoiStats.meanRad(q)=mean(R);        
                
    else %Set empty
        voronoiStats.ellipse(q,:)=nan;
        voronoiStats.ellipticity(q)=nan;
        voronoiStats.circularity(q)=nan;
        voronoiStats.angle(q)=nan;
        voronoiStats.area(q)=nan;
        voronoiStats.radii{q}=nan;
        voronoiStats.minRad(q)=nan;
        voronoiStats.maxRad(q)=nan;
        voronoiStats.meanRad(q)=nan;
    end
end

%% Compute histograms

E=voronoiStats.ellipticity;
E=E(~isnan(E));
hbE=linspace(min(E(:)),max(E(:)),numBins)';
hcE = histc(E,hbE);
hcE=hcE./size(E,1); %Normalize
voronoiStats.ellipticityHist.hc=hcE;
voronoiStats.ellipticityHist.hb=hbE;

E=voronoiStats.angle;
E=E(~isnan(E));
hbE=linspace(min(E(:)),max(E(:)),numBins)';
% hbE=linspace(0,pi,numBins)';
hcE = histc(E,hbE);
hcE=hcE./size(E,1); %Normalize
voronoiStats.angleHist.hc=hcE;
voronoiStats.angleHist.hb=hbE;

E=voronoiStats.area;
E=E(~isnan(E));
hbE=linspace(min(E(:)),max(E(:)),numBins)';
hcE = histc(E,hbE);
hcE=hcE./size(E,1); %Normalize
voronoiStats.areaHist.hc=hcE;
voronoiStats.areaHist.hb=hbE;

E=voronoiStats.circularity;
E=E(~isnan(E));
hbE=linspace(min(E(:)),max(E(:)),numBins)';
% hbE=linspace(0,1,numBins)';
hcE = histc(E,hbE);
hcE=hcE./size(E,1); %Normalize
voronoiStats.circularityHist.hc=hcE;
voronoiStats.circularityHist.hb=hbE;

E=voronoiStats.meanRad;
E=E(~isnan(E));
hbE=linspace(min(E(:)),max(E(:)),numBins)';
hcE = histc(E,hbE);
hcE=hcE./size(E,1); %Normalize
voronoiStats.meanRadHist.hc=hcE;
voronoiStats.meanRadHist.hb=hbE;

 
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
