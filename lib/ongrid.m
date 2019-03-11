function [L1t,L2t,L3t]=ongrid(V1,V2,V3,cm,vm)

L1t=ones(size(V1,1),size(cm,1));
L2t=ones(size(V2,1),size(cm,1));
L3t=ones(size(V3,1),size(cm,1));

for i=1:1:size(cm,1)
    
    c=cm(i,:);
    v=vm(i,:);
    
    [X1,Y1,Z1]=snap2grid(V1(:,1),V1(:,2),V1(:,3),c,v);
    X1=(X1-c(1))./v(1); Y1=(Y1-c(2))./v(2); Z1=(Z1-c(3))./v(3);
    V1s=[X1(:) Y1(:) Z1(:)];
    
    [X2,Y2,Z2]=snap2grid(V2(:,1),V2(:,2),V2(:,3),c,v);
    X2=(X2-c(1))./v(1); Y2=(Y2-c(2))./v(2); Z2=(Z2-c(3))./v(3);
    V2s=[X2(:) Y2(:) Z2(:)];
    
    [X3,Y3,Z3]=snap2grid(V3(:,1),V3(:,2),V3(:,3),c,v);
    X3=(X3-c(1))./v(1); Y3=(Y3-c(2))./v(2); Z3=(Z3-c(3))./v(3);
    V3s=[X3(:) Y3(:) Z3(:)];
    
    Vsnap=[V1s;V2s;V3s];
    Vsnap=1+(Vsnap-ones(size(Vsnap,1),1)*min(Vsnap,[],1));
    vsiz=max(Vsnap,[],1);
    Vsnap1=Vsnap(1:size(V1,1),:);
    Vsnap2=Vsnap(size(V1,1)+1:size(V1,1)+size(V2,1),:);
    Vsnap3=Vsnap(size(V1,1)+size(V2,1)+1:size(V1,1)+size(V2,1)+size(V3,1),:);
    
    INDsnap1=sub2ind(vsiz,Vsnap1(:,1),Vsnap1(:,2),Vsnap1(:,3));
    INDsnap2=sub2ind(vsiz,Vsnap2(:,1),Vsnap2(:,2),Vsnap2(:,3));
    INDsnap3=sub2ind(vsiz,Vsnap3(:,1),Vsnap3(:,2),Vsnap3(:,3));
    
    L1t(:,i)=ismembc(INDsnap1,sort(INDsnap2)) & ismembc(INDsnap1,sort(INDsnap3));
    L2t(:,i)=ismembc(INDsnap2,sort(INDsnap1)) & ismembc(INDsnap2,sort(INDsnap3));
    L3t(:,i)=ismembc(INDsnap3,sort(INDsnap1)) & ismembc(INDsnap3,sort(INDsnap2));
    
end

L1t=any(L1t,2); 
L2t=any(L2t,2); 
L3t=any(L3t,2);

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
