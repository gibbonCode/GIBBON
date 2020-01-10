function [V]=quad_smooth(V,Ev,IND_V_not_Ev,Lc,n,Ls,w1)

L1=Ev>0;
L2=IND_V_not_Ev>0;
w2=1-w1;

p=V;

for i=1:n
    Xp=NaN(size(Ev));  Yp=NaN(size(Ev));  Zp=NaN(size(Ev));
    Xp(L1)=p(Ev(L1),1); Yp(L1)=p(Ev(L1),2); Zp(L1)=p(Ev(L1),3);
    Vp1=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)];
    
    Xp=NaN(size(Ev));  Yp=NaN(size(Ev));  Zp=NaN(size(Ev));
    Xp(L2)=p(IND_V_not_Ev(L2),1); Yp(L2)=p(IND_V_not_Ev(L2),2); Zp(L2)=p(IND_V_not_Ev(L2),3);
    Vp2=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)];
    
    Vp=(w1.*Vp1+w2.*Vp2)./(w1+w2);
    p=p+Ls.*(Vp-p);
    
    p(Lc,:)=V(Lc,:);    
    if i>1
        D=gnansum((p-p_old).^2,2);
        SSQD_new=gnansum(D(:));
        if i>2
            SSQD_ratio=SSQD_new./SSQD_old;
            disp(num2str(SSQD_ratio));
        end
        SSQD_old=SSQD_new;
    end    
    p_old=p;
    
end
V=p;
 
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
