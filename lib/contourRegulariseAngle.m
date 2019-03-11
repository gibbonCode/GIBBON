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
        [~,ind1,ind2]=unique(pround(Vs,numDigitsMerge),'rows');
       
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
