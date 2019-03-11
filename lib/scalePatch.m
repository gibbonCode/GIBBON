function [Fc,Vc]=scalePatch(varargin)

% --------------------------------------------------------------------
% function [Fc,Vc]=scalePatch(F,V,scaleFactor)
%
% CHANGE LOG: 
% 19/12/2013 Fixed error related to single face entry, see if statement
% related to size(F,1)
% 2018/05/07 Created to be similar to patchDetach, which now calls this
% function
%
% --------------------------------------------------------------------
%%

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        scaleFactor=1;        
    case 3
        F=varargin{1};
        V=varargin{2};
        scaleFactor=varargin{3};
end

if numel(scaleFactor)==size(V,1) %If specified on the nodes
    scaleFactor=vertexToFaceMeasure(F,scaleFactor); % Convert to face metric
end

if numel(scaleFactor)~=1 && numel(scaleFactor)~=size(F,1)
    error('The number of elements in scaleFactor should be equal to 1 or the number of faces or the number of nodes');
end

%%
Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    if any(scaleFactor~=1)
        FX_mean=mean(FX,2);
        FX=((FX-FX_mean).*scaleFactor)+FX_mean;
    end
    Vc(:,q)=FX(:);
end
    
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));

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
