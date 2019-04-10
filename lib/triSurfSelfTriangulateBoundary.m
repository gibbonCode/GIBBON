function [F1,V1,ind1]=triSurfSelfTriangulateBoundary(F1,V1,ind1,angleThreshold,isClosedLoop)

% function [F1,V1,ind1]=triSurfSelfTriangulateBoundary(F1,V1,ind1,angleThreshold,isClosedLoop)
% ------------------------------------------------------------------------
%
%
%
% ------------------------------------------------------------------------

%%

ind1=ind1(:); %Force column

if ind1(1)==ind1(end)
    ind1=ind1(1:end-1);
    cropNeeded=1;
else
    cropNeeded=0;
end

%Get non-unique edge set
E=patchEdges(F1,0);

%Check if path order is consistent with edges
if ~any(all(E(:,1)==ind1(1) & E(:,2)==ind1(2),2))
    orderFlipped=1;
    ind1=flip(ind1); %Flip if not consistent    
else
    orderFlipped=0;
end

while 1
    
    [A]=patchPathAngles(F1,V1,ind1,isClosedLoop);
    
    indSharp=find(A<=angleThreshold);
    if isempty(indSharp)
        break
    end
    A_sharp=A(indSharp);
    
    if isClosedLoop==1
        es=[ind1(1:end-1) ind1(2:end); ind1(end) ind1(1); ]; %Closed loop
    else
        es=[ind1(1:end-1) ind1(2:end);]; %Open segment
    end
    
    [~,indNow]=min(A_sharp); %Sharpest first
    indVertexNow=ind1(indSharp(indNow));
    indEdges=(es(any(es==indVertexNow,2),:));
    f=fliplr(ind1(ismember(ind1,indEdges(:)))');

    if size(f,2)~=3
        warning('Invalid new triangle proposed');
        break
    end
    
    F1=[F1;f];    
    
    ind1=ind1(ind1~=indVertexNow);
    ind1=ind1(:);    
    
end

if orderFlipped
    ind1=flip(ind1);
end

if cropNeeded
    ind1(end+1)=ind1(1);
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
