function [Vi]=surfaceIntersect(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,interpMethod)

% function [Xi1,Yi1,Zi1]=surface_intersect(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
% ------------------------------------------------------------------------
% Determines the intersection points of 3 surfaces. X1, Y1, Z1 need to be a
% monotonic plaid surface.
%
% %%EXAMPLE
% clear all; close all; clc;
% 
% %% Simulating 3 surfaces
% n=75;
% xyzRange=linspace(-5,5,n);
% [X1,Y1]=meshgrid(xyzRange);
% Z1=(0.5.*sin(X1))+(0.5.*sin(Y1));
% [Z2,Y2]=meshgrid(xyzRange);
% X2=(0.5.*sin(Y2))+(0.5.*sin(Z2));
% [X3,Z3]=meshgrid(xyzRange);
% Yg=(1-(3.*exp(0.2.*(-(X3).^2-(Z3).^2))));
% Y3=(0.2.*sin(X3))+(0.2.*sin(Z3))+Yg;
% 
% %% Calculate intersection points
% 
% interpMethod='cubic';
% [V]=surfaceIntersect(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,interpMethod);
% 
% %% Plotting surfaces and intersection points
% figure; colordef(gcf,'white'); set(gcf,'outerposition',get(0,'ScreenSize'),'Color',[1 1 1]);  hold on; 
% ht1=xlabel('X'); ht2=ylabel('Y'); ht3=zlabel('Z'); ht4=title(['Found ',num2str(size(V,1)),' intersection point(s)']); set([ht1 ht2 ht3 ht4],'FontSize',20);
% plot3(V(:,1),V(:,2),V(:,3),'w.','MarkerSize',60);
% hs1=surf(X1,Y1,Z1,'FaceColor','r'); hs2=surf(X2,Y2,Z2,'FaceColor','g'); hs3=surf(X3,Y3,Z3,'FaceColor','b');
% set([hs1 hs2 hs3],'EdgeColor','none','FaceAlpha',0.7);view(3); axis equal; axis tight; axis vis3d; grid on; set(gca,'FontSize',20);
% camlight('headlight'); lighting phong; 
% drawnow;
% %% END EXAMPLE
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 30/08/2012
% ------------------------------------------------------------------------

%% Compute 2 sets of surface intersection curves
[Xi_12,Yi_12,Zi_12]=surfacePairIntersect(X1,Y1,Z1,X2,Y2,Z2,interpMethod); %Intersection curves of surface set 1 and 2
[Xi_13,Yi_13,~]=surfacePairIntersect(X1,Y1,Z1,X3,Y3,Z3,interpMethod); %Intersection curves of surface set 1 and 3

%% Computer intersection points for the two curve sets

Vi=[]; %Output array for unknown number of intersection points
for q1=1:numel(Xi_12)
    for q2=1:numel(Xi_13);
        Xi_12=Xi_12{q1}; Yi_12=Yi_12{q1}; Zi_12=Zi_12{q1};
        Xi_13=Xi_13{q2}; Yi_13=Yi_13{q2};                
        [xi,yi] = polyxpoly(Xi_12,Yi_12,Xi_13,Yi_13); %Finding intersection sites in XY projection of curves (this assumes that there are not multiple solutions in Z direction) 
        %         [xi,yi] = custom_polyxpoly(Xi_12,Yi_12,Xi_13,Yi_13); 
        if ~isempty(xi) %If intersection points are found
            [~,indUni,~] = unique(round(Xi_12*1e6)); %Quick fix to remove double entries (rounded to 6th decimal place)
            zi=interp1(Xi_12(indUni),Zi_12(indUni),xi,'cubic'); %interpolate to find z-coordinate, could base on 3D surface but here kept 1D for simplicity
            Vi(end+1,:)=[xi yi zi]; 
        end        
    end
end

end