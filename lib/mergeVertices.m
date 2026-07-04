function [varargout]=mergeVertices(varargin)

% function [F,V,ind1,ind2]=mergeVertices(F,V,numDigitsMerge)

%% Parse input

switch nargin    
    case 2
        F=varargin{1};
        V=varargin{2};
        numDigitsMerge=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        numDigitsMerge=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

cellMode=isa(F,'cell'); 

if isempty(numDigitsMerge)    
    D=patchEdgeLengths(F,V);
    if cellMode
        mean_D=mean(cellfun(@mean,D)); %Take mean across cell entries
    else
        mean_D=mean(D); %Mean across array
    end    
    numDigitsMerge=6-numOrder(mean_D); %base number of digits on mean
end

%% Merge nodes

[~,indKeep,indFix]=unique(pround(V,numDigitsMerge),'rows');
V=V(indKeep,:);

%% Fix indices in face array
if isa(F,'cell')
    for q=1:1:numel(F)
        F{q}=fixFaces(F{q},indFix);
    end
else
    F=fixFaces(F,indFix);
end

%% Collect output

varargout{1}=F;
varargout{2}=V;
varargout{3}=indKeep;
varargout{4}=indFix;

end

%% Fix indices in face array
function F=fixFaces(F,indFix)
    if size(F,1)==1
        F=indFix(F)'; %Fix indices in F
    else
        F=indFix(F); %Fix indices in F
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
