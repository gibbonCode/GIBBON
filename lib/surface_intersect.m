function [Vi]=surface_intersect(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,interpMethod)

% function [Xi1,Yi1,Zi1]=surface_intersect(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,interpMethod)
% ------------------------------------------------------------------------
% Determines the intersection points of 3 surfaces. X1, Y1, Z1 need to be a
% monotonic plaid surface.
%
% %%EXAMPLE
%
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 30/08/2012
% ------------------------------------------------------------------------

%% Compute intersection curves of surface set 1 and 2

[Xi_12,Yi_12,Zi_12]=surfacePairIntersect(X1,Y1,Z1,X2,Y2,Z2,interpMethod);

%% Compute intersection curves of surface set 1 and 3

[Xi_13,Yi_13,Zi_13]=surfacePairIntersect(X1,Y1,Z1,X3,Y3,Z3,interpMethod);

%% Computer intersection points for the two curve sets 12 and 13

Xi=[]; Yi=[]; Zi=[]; %Unknown number of intersection points
for i=1:numel(Xi_12)
    for j=1:numel(Xi_13);
        Xi_12=cell2mat(Xi_12(i,1));
        Yi_12=cell2mat(Yi_12(i,1));
        Zi_12=cell2mat(Zi_12(i,1));
        Xc2_i=cell2mat(Xi_13(j,1));
        Yc2_i=cell2mat(Yi_13(j,1));
        Zc2_i=cell2mat(Zi_13(j,1));

        %Finding intersection sites in XY projection of curves
        [xi,yi] = polyxpoly(Xi_12,Yi_12,Xc2_i,Yc2_i);
        Xi=[Xi; xi]; Yi=[Yi; yi];
        if ~isempty(xi)
            [Xi_12, m, n] = unique(Xi_12);
            Zi_12=Zi_12(m);
            zi=interp1(Xi_12,Zi_12,xi,'cubic');
            Zi=[Zi; zi];
        end
        
    end
end

Vi=[Xi Yi Zi];

end


%% OLD VERSION
% 
% 
% %% CHECK INPUT
% 
% % Check for plaid data.
% x1 = X1(1,:); y1 = Y1(:,1);
% if (size(X1,2)>1 && ~isequal(repmat(x1,size(X1,1),1),X1)) || ...
%         (size(Y1,1)>1 && ~isequal(repmat(y1,1,size(Y1,2)),Y1)),
%     error('MATLAB:interp2:meshgrid',...
%         ['X1 and Y1 must be matrices produced by MESHGRID. Use' ...
%         ' GRIDDATA instead \nof INTERP2 for scattered data.']);
% end
% 
% %% STEP 1:
% % L=~isnan(X1) & ~isnan(Y1) & ~isnan(Z1);
% % ZI_12 = griddata(X1(L),Y1(L),Z1(L),X2,Y2,'nearest');
% 
% ZI_12 = interp2(X1,Y1,Z1,X2,Y2,'linear');
% 
% C=(ZI_12-Z2);
% set(figure,'Visible','off');
% [D,h] = contour(X2,Y2,C,[0 0],'Visible','off','LineColor','none');
% 
% Xc1={};Yc1={};
% A1=get(h,'children');
% for i=1:1:length(A1)
%     B=get(A1(i));
%     Xc1_i=B.XData(~isnan(B.XData)); Xc1{i,1}=Xc1_i;
%     Yc1_i=B.YData(~isnan(B.YData)); Yc1{i,1}=Yc1_i;
%     Zc1_i=interp2(X1,Y1,Z1,Xc1_i,Yc1_i,'linear'); Zc1{i,1}=Zc1_i;
% %     plot3(Xc1_i(:), Yc1_i(:),Zc1_i(:),'k-','LineWidth',2);
% end
% 
% close gcf;
% %% STEP 2:
% 
% % ZI_13 = griddata(X1(L),Y1(L),Z1(L),X3,Y3,'nearest');
% ZI_13 = interp2(X1,Y1,Z1,X3,Y3,'linear');
% 
% C=(ZI_13-Z3);
% set(figure,'Visible','off');
% [D,h] = contour(X3,Y3,C,[0 0],'Visible','off','LineColor','none');
% 
% Xc2={};Yc2={};
% A2=get(h,'children');
% for i=1:1:length(A2)
%     B=get(A2(i));
%     Xc2_i=B.XData(~isnan(B.XData)); Xc2{i,1}=Xc2_i;
%     Yc2_i=B.YData(~isnan(B.YData)); Yc2{i,1}=Yc2_i;
%     Zc2_i=interp2(X1,Y1,Z1,Xc2_i,Yc2_i,'linear'); Zc2{i,1}=Zc2_i;
% %     plot3(Xc2_i(:), Yc2_i(:),Zc2_i(:),'k-','LineWidth',2);
% end
% 
% close gcf; 
% 
% %% STEP 3:
% 
% Xi1=[]; Yi1=[]; Zi1=[];
% Xi2=[]; Zi2=[];
% Yi3=[]; Zi3=[];
% for i=1:length(A1)
%     for j=1:length(A2);
%         Xc1_i=cell2mat(Xc1(i,1));
%         Yc1_i=cell2mat(Yc1(i,1));
%         Zc1_i=cell2mat(Zc1(i,1));
%         Xc2_i=cell2mat(Xc2(j,1));
%         Yc2_i=cell2mat(Yc2(j,1));
%         Zc2_i=cell2mat(Zc2(j,1));
% 
%         %Finding intersection sites in XY projection of curves
%         [xi,yi] = polyxpoly(Xc1_i,Yc1_i,Xc2_i,Yc2_i);
%         Xi1=[Xi1; xi]; Yi1=[Yi1; yi];
%         if ~isempty(xi)
%             [Xc1_i, m, n] = unique(Xc1_i);
%             Zc1_i=Zc1_i(m);
%             zi=interp1(Xc1_i,Zc1_i,xi,'cubic');
%             Zi1=[Zi1; zi];
%         end
%     end
% end

 
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
