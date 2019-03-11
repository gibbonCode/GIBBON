function [Xi,Yi,Zi]=surfacePairIntersect(X1,Y1,Z1,X2,Y2,Z2,interpMethod)

% function [Xi,Yi,Zi]=surfacePairIntersect(X1,Y1,Z1,X2,Y2,Z2,interpMethod)
% ------------------------------------------------------------------------
%
% This function employs interpolation and the CONTOUR command to compute 3D
% intersection curves for a set of two surfaces. Surfaces are required to
% be in matrix form i.e. monotonicly plaid. Surface 1 should ideally be of
% the form Z=f(X,Y), where f representes some suitable function. 
%
% Surface 1 is defined by X1, Y1 and Z1 
% Surface 2 is defined by X2, Y2 and Z2 
% interpMethod is the interpolation method used the following methods are
% supported: 'cubic', 'natural', 'linear' and 'nearest'. The latter three
% are based on scatteredInterpolant while the first is based on INTERP2. 
%
%
% %% EXAMPLE
% clear all; close all; clc;
% 
% %% Simulating complex surface pairs
% 
% %Surface 1
% n1=60;
% [X1,Y1]=meshgrid(linspace(-4,4,n1));
% Z1=X1+peaks(X1,Y1); Z1=pi.*Z1./max(Z1(:));
% 
% %Surface 2
% n2=75;
% [X2,Z2]=meshgrid(linspace(-pi,pi,n2));
% Y2=Z2+sin(2*X2)+sin(Z2);
% 
% %% Computer intersection curves
% [Xi,Yi,Zi]=surfacePairIntersect(X1,Y1,Z1,X2,Y2,Z2,'cubic');
% 
% %% Plotting surfaces and intersection curve
% 
% figure; colordef(gcf,'white'); set(gcf,'outerposition',get(0,'ScreenSize'),'Color',[1 1 1]);  hold on; 
% ht1=xlabel('X'); ht2=ylabel('Y'); ht3=zlabel('Z'); ht4=title(['Found ',num2str(numel(Xi)),' intersection curves']); set([ht1 ht2 ht3 ht4],'FontSize',20);
% h1=surf(X1,Y1,Z1,'FaceColor','r'); h2=surf(X2,Y2,Z2,'FaceColor','g');
% set([h1 h2],'EdgeColor','k','FaceAlpha',0.8);
% view(3); axis equal; axis tight; axis vis3d; grid on; set(gca,'FontSize',20);
% camlight('headlight'); lighting phong; 
% drawnow;
% pcolors=jet(numel(Xi)+1);
% for q=1:1:numel(Xi);
%     hp=plot3(Xi{q}, Yi{q},Zi{q},'k-','LineWidth',10); set(hp,'color',pcolors(q,:));
% end
% %% END EXAMPLE
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 30/08/2012
%------------------------------------------------------------------------

%Interpolate to find Z1-coordinates of surface 2 
if strcmp(interpMethod,'natural') || strcmp(interpMethod,'linear') || strcmp(interpMethod,'nearest') %TriScatterdInterp function
    F=scatteredInterpolant([X1(:) Y1(:)],Z1(:),interpMethod); %construct delaunay based interpolator
    Z2_I=F([X2(:) Y2(:)]); %Interpolate
    Z2_I=reshape(Z2_I,size(X2)); %This assumes no double points occured in delaunay tesselation
elseif strcmp(interpMethod,'cubic')
    % Data needs to be plaid and monotonically increasing
   Z2_I = interp2(X1,Y1,Z1,X2,Y2,interpMethod);
else
    error('myApp:argChk','Invalid interpolation method!');
end

%Calculate difference between two surfaces i.e. zero at intersection points
D=(Z2_I-Z2); 

%Derive the 2D "zero contour" i.e. the intersection contour to find the X
%and Y coordinates of the intersection curve(s)
set(figure,'Visible','off'); %open "invisible figure"
figureHandle=gcf;
[~,contourHandle] = contour(X2,Y2,D,[0 0],'Visible','off','LineColor','none'); %The "zero contours" 

%Get the contour object children (contains coordinates of contours)
contourChildren=get(contourHandle,'children'); 

%Get seperate intersection curves from contour 
Xi=cell(numel(contourChildren),1); Yi=Xi; Zi=Xi; %Cell arrays for output
for q=1:1:numel(contourChildren)   
    contourChild=get(contourChildren(q)); %Get child   
    xi=contourChild.XData(~isnan(contourChild.XData)); %X coordinate    
    yi=contourChild.YData(~isnan(contourChild.YData)); %Y coordinate
    
    %interpolate to obtain the Z-coordinate 
    if strcmp(interpMethod,'natural') || strcmp(interpMethod,'linear') || strcmp(interpMethod,'nearest') %TriScatterdInterp function
        F=scatteredInterpolant([X1(:) Y1(:)],Z1(:),interpMethod); %construct delaunay based interpolator
        zi=F([xi(:) yi(:)]); %Interpolate
        zi=reshape(zi,size(xi)); %This assumes no double points occured in delaunay tesselation
    else        
        zi=interp2(X1,Y1,Z1,xi,yi,interpMethod); 
    end
    %Store coordinates of current curve component in cell arrays
    Xi{q,1}=xi; Yi{q,1}=yi; Zi{q,1}=zi; 
end
close(figureHandle); %Close hidden figure

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
