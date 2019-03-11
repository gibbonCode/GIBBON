function [P]=trisurf_intersect(TRI1c,N1c,Vn1c,V1,TRI2c,N2c,Vn2c,V2,TRI3c,N3c,Vn3c,V3,d)


%% PREPARING INPUT 1
IND1n=repmat((1:numel(Vn1c(:,1)))'*ones(1,size(Vn2c,1)),[1 1 size(Vn3c,1)]);
IND1n=IND1n(:);

Ax=repmat(Vn1c(:,1)*ones(1,size(Vn2c,1)),[1 1 size(Vn3c,1)]);
Ay=repmat(Vn1c(:,2)*ones(1,size(Vn2c,1)),[1 1 size(Vn3c,1)]);
Az=repmat(Vn1c(:,3)*ones(1,size(Vn2c,1)),[1 1 size(Vn3c,1)]);
V1n=[Ax(:) Ay(:) Az(:)];

Ax=repmat(N1c(:,1)*ones(1,size(N2c,1)),[1 1 size(N3c,1)]);
Ay=repmat(N1c(:,2)*ones(1,size(N2c,1)),[1 1 size(N3c,1)]);
Az=repmat(N1c(:,3)*ones(1,size(N2c,1)),[1 1 size(N3c,1)]);
N1n=[Ax(:) Ay(:) Az(:)];

%% PREPARING INPUT 2
IND2n=repmat(ones(size(Vn1c,1),1)*(1:numel(Vn2c(:,1))),[1 1 size(Vn3c,1)]);
IND2n=IND2n(:);

Bx=repmat(ones(size(Vn1c,1),1)*Vn2c(:,1)',[1 1 size(Vn3c,1)]);
By=repmat(ones(size(Vn1c,1),1)*Vn2c(:,2)',[1 1 size(Vn3c,1)]);
Bz=repmat(ones(size(Vn1c,1),1)*Vn2c(:,3)',[1 1 size(Vn3c,1)]);
V2n=[Bx(:) By(:) Bz(:)];

Bx=repmat(ones(size(N1c,1),1)*N2c(:,1)',[1 1 size(N3c,1)]);
By=repmat(ones(size(N1c,1),1)*N2c(:,2)',[1 1 size(N3c,1)]);
Bz=repmat(ones(size(N1c,1),1)*N2c(:,3)',[1 1 size(N3c,1)]);
N2n=[Bx(:) By(:) Bz(:)];

%% PREPARING INPUT 3
IND3n=permute(repmat((1:numel(Vn3c(:,1)))'*ones(1,size(Vn1c,1)),[1 1 size(Vn2c,1)]),[2 3 1]);
IND3n=IND3n(:);

Cx=permute(repmat(Vn3c(:,1)*ones(1,size(Vn1c,1)),[1 1 size(Vn2c,1)]),[2 3 1]);
Cy=permute(repmat(Vn3c(:,2)*ones(1,size(Vn1c,1)),[1 1 size(Vn2c,1)]),[2 3 1]);
Cz=permute(repmat(Vn3c(:,3)*ones(1,size(Vn1c,1)),[1 1 size(Vn2c,1)]),[2 3 1]);
V3n=[Cx(:) Cy(:) Cz(:)];

Cx=permute(repmat(N3c(:,1)*ones(1,size(N1c,1)),[1 1 size(N2c,1)]),[2 3 1]);
Cy=permute(repmat(N3c(:,2)*ones(1,size(N1c,1)),[1 1 size(N2c,1)]),[2 3 1]);
Cz=permute(repmat(N3c(:,3)*ones(1,size(N1c,1)),[1 1 size(N2c,1)]),[2 3 1]);
N3n=[Cx(:) Cy(:) Cz(:)];

%% Find triangle plane intersections
X=plane_intersect(V1n,V2n,V3n,N1n,N2n,N3n);
L_notnan=~any(isnan(X),2);

%% Check if intersection points lie on triangle
if isempty(X) || all(~L_notnan)
    P=[];
else    
    X=X(L_notnan,:);    
    L1=intriangle(TRI1c(IND1n(L_notnan),:),N1c(IND1n(L_notnan),:),V1,X,d);
    L2=intriangle(TRI2c(IND2n(L_notnan),:),N2c(IND2n(L_notnan),:),V2,X,d);
    L3=intriangle(TRI3c(IND3n(L_notnan),:),N3c(IND3n(L_notnan),:),V3,X,d);
    L=L1&L2&L3;    
    P=X(L,:);
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
