function [varargout]=meshBoundary(varargin)

switch nargin
    case 1
        E=varargin{1};
        [F,CF]=element2patch(E,[]); %Get mesh faces
    case 2
        E=varargin{1};
        elementType=varargin{2};
        [F,CF]=element2patch(E,[],elementType); %Get mesh faces
    case 3        
        E=varargin{1};        
        elementType=varargin{2};        
        C=varargin{3};
        [F,CF]=element2patch(E,C,elementType); %Get mesh faces
    otherwise
        error('Wrong number of input arguments');
end

if isempty(CF)
    CF=zeros(size(F,1),1);
end

%Check shared faces faces
numFacesIni=size(F,1);
[F_uni,indF,IND_F_2]=uniqueIntegerRow(F);

CF=CF(indF,:);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Compose boundary set from faces that are used once
logicBoundary=F_count==1; 
Fb=F_uni(logicBoundary,:);
CFb=CF(logicBoundary,:);

varargout{1}=Fb;
varargout{2}=F_count;
varargout{3}=CFb;
 
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
