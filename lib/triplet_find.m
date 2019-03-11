function [L1t,L2t,L3t]=triplet_find(DT1,DT2,DT3,Vn1,Vn2,Vn3,T_dist,opt)

switch opt
    case 1 %Full Delaunay input based
        [~,Dc12] = nearestNeighbor(DT2,Vn1);
        [~,Dc13] = nearestNeighbor(DT3,Vn1);
        L1t=Dc12<=T_dist & Dc13<=T_dist;
        
        [~,Dc21] = nearestNeighbor(DT1,Vn2);
        [~,Dc23] = nearestNeighbor(DT3,Vn2);
        L2t=Dc21<=T_dist & Dc23<=T_dist;
        
        [~,Dc31] = nearestNeighbor(DT1,Vn3);
        [~,Dc32] = nearestNeighbor(DT2,Vn3);
        L3t=Dc31<=T_dist & Dc32<=T_dist;
        
    case 2 %Accelerating renew with cropped Delaunay if possible otherwise full
        
        [~,Dc12] = nearestNeighbor(DT2,Vn1); 
        L12t=Dc12<=T_dist;
        if any(L12t)
            [~,Dc13] = nearestNeighbor(DT3,Vn1);
            L1t=L12t & Dc13<=T_dist;
        else
            L1t=false(size(Vn1,1),1);
        end
        
        if any(L1t)            
            if nnz(L1t)>=4
                DT=delaunayTriangulation(Vn1(L1t,:)); %Tesselation of selected coordinates
                if ~isempty(DT.ConnectivityList)
                    DT1=DT;
                end
            end
            
            [~,Dc21] = nearestNeighbor(DT1,Vn2);
            L21t=Dc21<=T_dist;
            if any(L21t)
                [~,Dc23] = nearestNeighbor(DT3,Vn2);            
                L2t=L21t & Dc23<=T_dist;
            else
                L2t=false(size(Vn2,1),1);
            end
        else
            L2t=false(size(Vn2,1),1);
        end
        
        if any(L2t)
            if nnz(L2t)>=4
                DT=delaunayTriangulation(Vn2(L2t,:)); %Tesselation of selected coordinates
                if ~isempty(DT.ConnectivityList)
                    DT2=DT;
                end
            end
            
            [~,Dc31] = nearestNeighbor(DT1,Vn3);
            L31t=Dc31<=T_dist;
            if any(L31t)
                [~,Dc32] = nearestNeighbor(DT2,Vn3);
                L3t=L31t & Dc32<=T_dist;
            else
                L3t=false(size(Vn3,1),1);
            end
        else
            L3t=false(size(Vn3,1),1);
        end
        
    case 3 %Accelerating renew with cropped Delaunay/non-Delaunay
        [~,Dc12] = nearestNeighbor(DT2,Vn1);
        [~,Dc13] = nearestNeighbor(DT3,Vn1);
        L1t=Dc12<=T_dist & Dc13<=T_dist;
        if any(L1t)
            if nnz(L1t)>=4
                DT=delaunayTriangulation(Vn1(L1t,:)); %Tesselation of selected coordinates
                if ~isempty(DT.ConnectivityList)
                    DT1=DT;
                end
                [~,Dc21] = nearestNeighbor(DT1,Vn2);
                [~,Dc31] = nearestNeighbor(DT1,Vn3);
            else
                [INDc,Dc12] = dsearchn(Vn2,Vn1(L1t,:));
                Dc21=NaN(size(Vn2,1),1); Dc21(INDc)=Dc12;
                
                [INDc,Dc13] = dsearchn(Vn3,Vn1(L1t,:));
                Dc31=NaN(size(Vn3,1),1); Dc31(INDc)=Dc13;
            end
            [~,Dc23] = nearestNeighbor(DT3,Vn2);
            L2t=Dc21<=T_dist & Dc23<=T_dist;
        else
            L2t=false(size(Vn2,1),1);
        end
        if any(L2t)
            if nnz(L2t)>=4 && nnz(L1t)>=4
                 DT=delaunayTriangulation(Vn2(L2t,:)); %Tesselation of selected coordinates
                if ~isempty(DT.ConnectivityList)
                    DT2=DT;
                end
                [~,Dc32] = nearestNeighbor(DT2,Vn3);
            else
                [INDc,Dc23] = dsearchn(Vn3,Vn2(L2t,:));
                Dc32=NaN(size(Vn3,1),1); Dc32(INDc)=Dc23;
            end
            L3t=Dc31<=T_dist & Dc32<=T_dist;
        else
            L3t=false(size(Vn3,1),1);
        end
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
