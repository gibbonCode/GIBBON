function [a,d]=vectorOrthogonalPair(f)

% function [a,d]=vectorOrthogonalPair(f)
% ------------------------------------------------------------------------
%
% Based on the input vector f this function generates the normalized output
% vectors a and d which are orthogonal to eachother and to f. 
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2017/05/08: Updated for GIBBON
% 2017/05/08: Fixed bug in relation to co-linear output.

%%

f=vecnormalize(f);

e1=[1 0 0];
e2=[0 1 0];
e3=[0 0 1];

E1=e1(ones(size(f,1),1),:);
E2=e2(ones(size(f,1),1),:);
E3=e3(ones(size(f,1),1),:);
 
% T=abs(real([acos(dot(f,E1,2)) acos(dot(f,E2,2)) acos(dot(f,E3,2))]));
% T(T>pi)=(T(T>pi)-pi);
% T=abs(T-0.5*pi);
% [~,J_min]=min(T,[],2);

[~,J_min]=min(abs(f),[],2);

d=zeros(size(f));
for q=1:1:3
   logicNow=J_min==q;   
   if any(logicNow)
       switch q
           case 1
               E=E1(logicNow,:);
           case 2
               E=E2(logicNow,:);
           case 3
               E=E3(logicNow,:);
       end
       d(logicNow,:)=cross(E,f(logicNow,:),2);
   end
end
[d]=vecnormalize(d);
a=cross(f,d,2); [a]=vecnormalize(a); %d is orthogonal to f and a
d=cross(f,a,2); [d]=vecnormalize(d); %a is reset to be orthogonal to both f and d
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
