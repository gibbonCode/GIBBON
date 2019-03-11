function [X,Y,Z]=contourRegularise(varargin)

switch nargin
    case 2
        Vcs=varargin{1};
        np=varargin{2};
        interpMethod='pchip'; %Default method
    case 3
        Vcs=varargin{1};
        np=varargin{2};
        interpMethod=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

%%
numContours=numel(Vcs);

%Allocate coordinate matrices
X=nan(numContours,np); Y=nan(numContours,np); Z=nan(numContours,np);
startDefined=0;
for q=1:1:numContours
    
    %Join goups if present N.B. assumes they belong to the same curve!
    Vs=[];
    for qGroup=1:1:numel(Vcs{q})
        Vss=Vcs{q}{qGroup};
        if ~isempty(Vss)
            if isPolyClockwise(Vss)
                Vss=flipud(Vss);
            end
            Vs=[Vs; Vss];
        end
    end
    
    %Resample curve so they have the same number of points
    if ~isempty(Vs)
        numDigitsMerge=6-numOrder(mean(sum(diff(Vs,[],1).^2,2)));
        [~,ind1,~]=unique(pround(Vs,numDigitsMerge),'rows');
        
        Vs=Vs(ismember(1:size(Vs,1),ind1),:); %This maintains point order
        
        [Vs]=evenlySampleCurve(Vs,np,interpMethod,1);
        
        if startDefined==0 %After 1 contour has been added
            V_start=Vs(1,:);
            startDefined=1;
        else
            %Reorder curves so start and end points match closely
            D=sqrt(sum((Vs-V_start(ones(1,size(Vs,1),1),:)).^2,2));
            [~,indMin]=min(D);
            if indMin~=1
                Vs=[Vs(indMin:size(Vs,1),:); Vs(1:indMin-1,:)];
            end
            V_start=Vs(1,:);
        end
        
        %Store coordinates
        [X(q,:),Y(q,:),Z(q,:)]=getColumns(Vs);
    end
end

%%
%Remove twist
for q=2:1:size(numContours,1)
    v1=[X(q-1,:)' Y(q-1,:)' Z(q-1,:)'];
    v2=[X(q,:)' Y(q,:)' Z(q,:)'];
    [v2f,~,~]=minPolyTwist(v1,v2);
    X(q,:)=v2f(:,1);
    Y(q,:)=v2f(:,2);
    Z(q,:)=v2f(:,3);
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
