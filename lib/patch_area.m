function [A]=patch_area(F,V)

% function [A]=patch_area(F,V)
% ------------------------------------------------------------------------
%This simple function calculates the areas of the faces specified by F and
%V. The output is a vector A containing size(F,1) elements. The face areas
%are calculated via triangulation of the faces. If faces are already
%triangular triangulation is skipped are area calculation is direction
%performed.
%
%%% EXAMPLE
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 12/04/2011
%------------------------------------------------------------------------
%%

if size(V,2)~=3
    V(:,3)=0;
end

if size(F,2)>3 %Convert to triangles    
    %Format of column index in F
    EColumnInd=[(1:size(F,2)); (1:size(F,2))];
    EColumnInd=[EColumnInd(2:end) EColumnInd(1)];
    
    %Derive edges matrix
    E=F(:,EColumnInd)'; %Use index into F to create edges matrix
    E=reshape(E,2,numel(E)/2)';
    
    %Add centre point to vertices
    X=V(:,1); Y=V(:,2); Z=V(:,3);
    if size(F,1)==1 %If its only a single face the behaviour of indexing changes
        Vc=[mean(X(F(:))) mean(Y(F(:))) mean(Z(F(:)))];
    else
        Vc=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];
    end
    Vt=[V;Vc];
    
    Ft=E;
    ind_face=ones(size(F,2),1)*(1:1:size(F,1));
    ind_face=ind_face(:);
    ind_centres=ind_face+size(V,1);
    Ft(:,3)=ind_centres;
    
    %Calculate surface areas
    V12=Vt(Ft(:,2),:)-Vt(Ft(:,1),:);
    V13=Vt(Ft(:,3),:)-Vt(Ft(:,1),:);
    A=abs(0.5.*sqrt(sum(cross(V12,V13,2).^2,2)));
    A=sum(reshape(A,size(F,2),size(F,1)),1)';
    
else %Already triangles
    %Calculate surface areas
    V12=V(F(:,2),:)-V(F(:,1),:);
    V13=V(F(:,3),:)-V(F(:,1),:);
    A=abs(0.5.*sqrt(sum(cross(V12,V13,2).^2,2)));
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
