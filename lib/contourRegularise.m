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
for q=1:1:numContours;
    
    %Join goups if present N.B. assumes they belong to the same curve!
    Vs=[];
    for qGroup=1:1:numel(Vcs{q})
        Vss=Vcs{q}{qGroup};
        if ~isempty(Vss)
            if isPolyClockwise(Vss);
                Vss=flipud(Vss);
            end
            Vs=[Vs; Vss];
        end
    end
    
    %Resample curve so they have the same number of points
    if ~isempty(Vs)
        
        [~,ind1,~]=unique(pround(Vs,5),'rows');
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
for q=2:1:size(numContours,1);
    v1=[X(q-1,:)' Y(q-1,:)' Z(q-1,:)'];
    v2=[X(q,:)' Y(q,:)' Z(q,:)'];
    [v2f,~,~]=minPolyTwist(v1,v2);
    X(q,:)=v2f(:,1);
    Y(q,:)=v2f(:,2);
    Z(q,:)=v2f(:,3);
end

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
