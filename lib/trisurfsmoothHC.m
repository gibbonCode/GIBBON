function [Vs]=trisurfsmoothHC(F,V,IND_V,alp,bet,tol,n,disp_on)

%Humphreys Classes Smoothening

% alpha [0..1] influence of original/previous points 0-> full previous, 1-> full original
% beta  [0..1] e.g. > 0.5

if isempty(IND_V)
    [~,IND_V]=patchIND(F,V);
end

L=IND_V>0;
Xp=NaN(size(IND_V));  Yp=NaN(size(IND_V));  Zp=NaN(size(IND_V));

p=V; cc=0; i=1;
while cc==0
    q=p;
    Xp(L)=p(IND_V(L),1); Yp(L)=p(IND_V(L),2); Zp(L)=p(IND_V(L),3);
    p=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of neighbourhood, Laplacian operation
    SSQD_new=gnansum((p(:)-q(:)).^2);
    if i>1
        SSQD_ratio=SSQD_new./SSQD_old;
        if disp_on
            disp(['Iteration ',num2str(i),', SSQD ratio: ',num2str(SSQD_ratio)]);
        end
        if abs(1-SSQD_ratio)<=tol
            cc=1; %Convergence tolerance achieved
            if disp_on
                disp('Convergence tolerance on SSQD ratio reached!');
            end
        end
    else
        if disp_on
            disp(['Iteration ',num2str(i),', SSQD initial: ',num2str(SSQD_new)]);
        end
    end
    
    if i>=n
        cc=1;
        if disp_on
            disp('Maximum number of iteration reached!');
        end
    end
    SSQD_old=SSQD_new;
    
    b=p-((alp.*V)+((1-alp).*q)); %Difference at centre vertex
    Xp(L)=b(IND_V(L),1); Yp(L)=b(IND_V(L),2); Zp(L)=b(IND_V(L),3);
    c=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of difference at centre vertex of neighbourhood
    p=p-((bet.*b)+(1-bet).*(c));
    i=i+1;
end

Vs=p;

 
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
