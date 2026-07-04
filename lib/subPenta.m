function [varargout]=subPenta(varargin)

% function [E,V,C,CV]=subPenta(E,V,n,splitMethod)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
%
% 2020/01/29 Created
% ------------------------------------------------------------------------


%% Parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        n=1;
        splitMethod=1;
    case 3
        E=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=1;
    case 4
        E=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=varargin{4};
end

C=(1:1:size(E,1))'; %Element colors or indices
CV=[];

%%
if n>0
    for qIter=1:1:n
        switch splitMethod
            case 1
                [E,V,C]=subPenta(E,V,n,2);
                [Es,Vs]=subPenta(E,V,n,3);                                
            case 2
                F1=E(:,[1 2 3]);
                F2=E(:,[4 5 6]);
                [F1s,V1s]=subtri(F1,V);
                [F2s,V2s]=subtri(F2,V);
                Vs=[V1s;V2s];
                Es=[F1s F2s+size(V1s,1)];                
            case 3
                %% Mid edge sets
                edgeMat=[E(:,[1 4]); E(:,[2 5]);  E(:,[3 6]);];   %top-bottom edges
                
                E_sort=sort(edgeMat,2); %Sorted edges matrix
                [~,ind1,~]=unique(E_sort,'rows');
                edgeMat=edgeMat(ind1,:);
                
                numPoints = size(V,1);
                numEdges = size(edgeMat,1);
                
                % Get indices of the four edges associated with each face
                A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
                A = max(A,A'); %Copy symmetric
                
                %Indices for A matrix
                indA_14=E(:,1)+(E(:,4)-1)*numPoints;
                indA_25=E(:,2)+(E(:,5)-1)*numPoints;
                indA_36=E(:,3)+(E(:,6)-1)*numPoints;
                
                %Get indices for vertex array
                indV_14=full(A(indA_14));
                indV_25=full(A(indA_25));
                indV_36=full(A(indA_36));
                
                %% Create element array
                Es=[E(:,1:3) indV_14 indV_25 indV_36;...
                    indV_14 indV_25 indV_36 E(:,4:6)];
                
                %% Create vertex array
                Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
                Vs = [V; Vn;]; %Join point sets
                CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1);];
        end
        
        %% Override input
        C=repmat(C,[size(Es,1)/size(E,1),1]);
        E=Es;
        V=Vs;
    end
end

varargout{1}=E;
varargout{2}=V;
varargout{3}=C;
varargout{4}=CV;


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
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
