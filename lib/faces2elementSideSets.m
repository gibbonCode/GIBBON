function [varargout]=faces2elementSideSets(varargin)
% function [elementId,elementSides,elemSets,elemSetSides]=faces2elementSideSets(f,E)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% 2020/05/30 Created
% ------------------------------------------------------------------------

%%

switch nargin
    case 2
        f=varargin{1};
        E=varargin{2};
        reorderOpt=0;
    case 3
        f=varargin{1};
        E=varargin{2};
        reorderOpt=varargin{3};
end
%%
logicTouch=sum(ismember(E,f),2)==size(f,2);
ind=(1:1:size(E,1))';
[Fs,C,CF]=element2patch(E(logicTouch,:),ind(logicTouch));
logicMember=isrowmember(Fs,f,reorderOpt);
elementId=C(logicMember);
elementSides=CF(logicMember);
varargout{1}=elementId;
varargout{2}=elementSides;
if nargout>2
    setSideUni=unique(elementSides);
    numSets=numel(setSideUni);
    elemSets=cell(1,numSets);
    elemSetSides=zeros(1,numSets);
    for q=1:1:numSets
        sideType=setSideUni(q);
        logicNow=elementSides==sideType;
        indElem=unique(elementId(logicNow));
        elemSets{q}=indElem(:)';
        elemSetSides(q)=sideType;
    end
    varargout{3}=elemSets;
    varargout{4}=elemSetSides;
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
