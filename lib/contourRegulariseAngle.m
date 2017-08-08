function [Xcs,Ycs,Zcs]=contourRegulariseAngle(varargin)

switch nargin
    case 2
        Vcs=varargin{1};
        np=varargin{2};
        thetaStart=0; %Default angle
        interpMethod='pchip'; %Default method
    case 3
        Vcs=varargin{1};
        np=varargin{2};
        thetaStart=varargin{3}; 
        interpMethod='pchip'; %Default method
            case 4
        Vcs=varargin{1};
        np=varargin{2};
        thetaStart=varargin{3}; 
        interpMethod=varargin{4}; %Default method
    otherwise        
        error('Wrong number of input arguments');
end

%% CONTROL PARAMETERS
closeLoopOpt=1;

numContours=numel(Vcs);

%Allocate coordinate matrices
Xcs=nan(numContours,np); Ycs=nan(numContours,np); Zcs=nan(numContours,np);
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
 
        [~,ind1,ind2]=unique(pround(Vs,5),'rows');
        Vs=Vs(ismember(1:size(Vs,1),ind1),:);
        
%         figure; plotV(Vs,'k.-'); hold on; 
%         plotV(Vs(1,:),'r.');
%         plotV(Vs(end,:),'b.-');
        
        [Vs]=evenlySampleCurve(Vs,np,interpMethod,closeLoopOpt);
        
        %Reorder curve so that start points are close
         %%
%         if startDefined==0 %After 1 contour has been added
%             [thetaAll,~,~] = cart2pol(Vs(:,1),Vs(:,2),Vs(:,3));
%             thetaDiff=thetaAll-thetaStart;
%             [~,indMin]=min(thetaDiff);
%             V_start=Vs(indMin,:);
%             %Reorder first curve
%             Vs=[Vs(indMin:size(Vs,1),:); Vs(1:indMin-1,:)];
%             startDefined=1;
%         else
%             %Reorder curves so start and end points match closely
%             D=sqrt(sum((Vs-V_start(ones(1,size(Vs,1),1),:)).^2,2));
%             [~,indMin]=min(D);
%             if indMin~=1
%                 Vs=[Vs(indMin:size(Vs,1),:); Vs(1:indMin-1,:)];
%             end
%             V_start=Vs(1,:);
%         end
        
        %%
        meanVs=mean(Vs,1);
        Vsc=Vs-meanVs(ones(size(Vs,1),1),:);
        [thetaAll,~,~] = cart2pol(Vsc(:,1),Vsc(:,2),Vsc(:,3));
        thetaDiff=thetaAll-thetaStart;
        [~,indMin]=min(thetaDiff);
        
        if indMin~=1
            Vs=[Vs(indMin:size(Vs,1),:); Vs(1:indMin-1,:)];
        end
        
        %%
        
        %Store coordinates
        [Xcs(q,:),Ycs(q,:),Zcs(q,:)]=getColumns(Vs);
%         any(isnan(Zcs(q,:)))
    end
end
 
%% 
% ********** _license boilerplate_ **********
% 
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
