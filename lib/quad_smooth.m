function [V]=quad_smooth(V,Ev,IND_V_not_Ev,Lc,n,Ls,w1)

L1=Ev>0;
L2=IND_V_not_Ev>0;
w2=1-w1;

p=V;

for i=1:n;
    Xp=NaN(size(Ev));  Yp=NaN(size(Ev));  Zp=NaN(size(Ev));
    Xp(L1)=p(Ev(L1),1); Yp(L1)=p(Ev(L1),2); Zp(L1)=p(Ev(L1),3);
    Vp1=[nanmean(Xp,2) nanmean(Yp,2) nanmean(Zp,2)];
    
    Xp=NaN(size(Ev));  Yp=NaN(size(Ev));  Zp=NaN(size(Ev));
    Xp(L2)=p(IND_V_not_Ev(L2),1); Yp(L2)=p(IND_V_not_Ev(L2),2); Zp(L2)=p(IND_V_not_Ev(L2),3);
    Vp2=[nanmean(Xp,2) nanmean(Yp,2) nanmean(Zp,2)];
    
    Vp=(w1.*Vp1+w2.*Vp2)./(w1+w2);
    p=p+Ls.*(Vp-p);
    
    p(Lc,:)=V(Lc,:);    
    if i>1
        D=nansum((p-p_old).^2,2);
        SSQD_new=nansum(D(:));
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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
