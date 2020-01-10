function [Vs]=trisurfsmooth(F,V,n,IND_V,par,opt)

% alpha [0..1] influence of original/previous points 0-> full previous, 1-> full original
% beta  [0..1] e.g. > 0.5 

if isempty(IND_V)
    [~,IND_V]=patchIND(F,V);
end

L=IND_V>0;
Xp=NaN(size(IND_V));  Yp=NaN(size(IND_V));  Zp=NaN(size(IND_V));

switch opt
    case 1 %Laplacian
        p=V;      
        Ls=par;
        for i=1:n         
            Xp(L)=p(IND_V(L),1); Yp(L)=p(IND_V(L),2); Zp(L)=p(IND_V(L),3);
            Vp=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)];
            p=p+Ls.*(Vp-p);    
        end
    case 2 %Humphreys Classes Smoothening
        p=V;        
        alp=par(1);
        bet=par(2);
        for i=1:n            
            q=p;
            Xp(L)=p(IND_V(L),1); Yp(L)=p(IND_V(L),2); Zp(L)=p(IND_V(L),3);            
            p=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of neighbourhood, Laplacian operation   
            SSQD_new=sum((p(:)-q(:)).^2);
            if i>1
                cf=SSQD_new./SSQD_old; disp(['SSQD ratio: ',num2str(cf)]);
            end
            SSQD_old=SSQD_new;
            
            b=p-((alp.*V)+((1-alp).*q)); %Difference at centre vertex
            
            Xp(L)=b(IND_V(L),1); Yp(L)=b(IND_V(L),2); Zp(L)=b(IND_V(L),3);            
            c=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of difference at centre vertex of neighbourhood           
            p=p-((bet.*b)+(1-bet).*(c));            
        end        
    case 3 %Humphreys Classes Smoothening convergence, n(1) is max number of iterations, n(2) is used as tolerance on SSQD ratio
        p=V;
        alp=par(1);
        bet=par(2);
        cc=0; i=1;
        while cc==0
            q=p;
            Xp(L)=p(IND_V(L),1); Yp(L)=p(IND_V(L),2); Zp(L)=p(IND_V(L),3);            
            p=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of neighbourhood, Laplacian operation   
            SSQD_new=sum((p(:)-q(:)).^2);
            if i>1
                SSQD_ratio=SSQD_new./SSQD_old; disp(['Iteration ',num2str(i),', SSQD ratio: ',num2str(SSQD_ratio)]);
                if abs(1-SSQD_ratio)<=n(2)
                    cc=1; %Convergence tolerance achieved
                    disp('Convergence tolerance on SSQD ratio reached!');
                end
            else
                disp(['Iteration ',num2str(i),', SSQD initial: ',num2str(SSQD_new)]);
            end
            
            if i>=n(1)
                cc=1;
                disp('Maximum number of iteration reacheds!');
            end
            SSQD_old=SSQD_new;            
            
            b=p-((alp.*V)+((1-alp).*q)); %Difference at centre vertex            
            Xp(L)=b(IND_V(L),1); Yp(L)=b(IND_V(L),2); Zp(L)=b(IND_V(L),3);            
            c=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)]; %Mean of difference at centre vertex of neighbourhood           
            p=p-((bet.*b)+(1-bet).*(c));
            i=i+1;
        end        
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
