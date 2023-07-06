function [Fn,Vn]=subTriSplit(varargin)

% function [Fn,Vn]=subTriSplit(F,V,n)
% ------------------------------------------------------------------------
% This function splits the triangulation defined by the faces F, and the
% vertices V, n times. Each triangle is linearly split into 4 triangles
% with each iterations. 
% 
% 2023/06/06 KMM: Added to GIBBON
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        n=1;        
    case 3
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
end

%% Process sub-triangulation through splitting

if n==0 % No action needed if n is 0
    Fn=F;
    Vn=V;
elseif n==1 % Splitting just once
    %Setting up new vertices and faces such that no new unshared points are
    % added. 
    
    E=[F(:,[1 2]); F(:,[2 3]);  F(:,[3 1])]; %Non-unique edges matrix
    Es=sort(E,2); %Sorted edges matrix, so 1-2 is the same as 2-1
    [~,ind1,ind2]=unique(Es,'rows'); %Get indices of unique edges
    E=E(ind1,:); %Unique edges        
    ind2=ind2+size(V,1); %Offset indices since new pointsd (Vm below) are appended "under" V

    %Create face array (knowing each face has three new points associated
    %with it)
    numFaces=size(F,1); % The number of faces
    indF  = 1:1:numFaces; % Indices for all faces
    indP1 = ind2(indF); % Indices of all new first points
    indP2 = ind2(indF+numFaces); % Indices of all new second points
    indP3 = ind2(indF+2*numFaces); % Indices of all new third points
        
    Fn = [F(:,1) indP1 indP3;... % 1st new corner face
          F(:,2) indP2 indP1;... % 2nd new corner face
          F(:,3) indP3 indP2;... % 3rd new corner face
          indP1 indP2 indP3];    % New central face

    %Create vertex arrays
    Vm=0.5*(V(E(:,1),:)+V(E(:,2),:)); %new mid-edge points
    Vn = [V; Vm]; %Join point sets    
elseif n>1 %Splitting more than once so recursively split once 
    Fn=F; Vn=V; %Initialise Fn, Vn as input F and V
    for q=1:1:n %Split once n times
        [Fn,Vn]=subTriSplit(Fn,Vn,1); %Split once
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