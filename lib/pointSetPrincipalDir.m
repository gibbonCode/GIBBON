function [varargout]=pointSetPrincipalDir(X)

% function [V,S,U]=pointSetPrincipalDir(X)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%%

%Cope with 2D input
if size(X,2)==2
    X(:,3)=0; %Force 3D
end

%Centre on own mean
MU=mean(X,1); %Point set mean
X=X-MU(ones(size(X,1),1),:); %Centre points around mean

%Compute singular value decomposition to get principal directions
[U,S,V]=svd(X,0); 

%% Collect output
varargout{1}=V;
varargout{2}=S;
varargout{3}=U;
 
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
