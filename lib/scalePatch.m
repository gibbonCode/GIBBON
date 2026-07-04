function [varargout]=scalePatch(varargin)

% ------------------------------------------------------------------------
% function [Fc,Vc]=scalePatch(F,V,scaleFactor,V_scale)
%
%
% CHANGE LOG: 
% 19/12/2013 Fixed error related to single face entry, see if statement
% related to size(F,1)
% 2018/05/07 Created to be similar to patchDetach, which now calls this
% function
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        scaleFactor=1;
        V_scale=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        scaleFactor=varargin{3};
        V_scale=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        scaleFactor=varargin{3};
        V_scale=varargin{4};
end

if numel(scaleFactor)==size(V,1) %If specified on the nodes
    scaleFactor=vertexToFaceMeasure(F,scaleFactor); % Convert to face metric
end

if numel(scaleFactor)~=1 && numel(scaleFactor)~=size(F,1)
    error('The number of elements in scaleFactor should be equal to 1 or the number of faces or the number of nodes');
end

%%
Vc=zeros(size(F,1)*size(F,2),size(V,2));
if nargout==3
    Vcc=Vc;
end

for q=1:1:size(V,2)
    X=V(:,q);
    
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    
    if nargout==3
        Vcc(:,q)=FX(:);
    end
    
    if any(scaleFactor~=1)
        if isempty(V_scale)
            FX_mean=mean(FX,2);
        else
            FX_mean=V_scale(:,q);
        end
        FX=((FX-FX_mean(:,ones(size(FX,2),1))).*scaleFactor)+FX_mean(:,ones(size(FX,2),1));
    end
    Vc(:,q)=FX(:);
end
    
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));

%% Collect output
varargout{1}=Fc;
varargout{2}=Vc;
if nargout==3
    varargout{3}=Vcc;
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
