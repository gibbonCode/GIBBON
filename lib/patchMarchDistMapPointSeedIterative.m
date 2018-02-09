function [indSeed,D_map,seedIndex]=patchMarchDistMapPointSeedIterative(F,V,indStart,numSeeds,W,optStruct)

%%

if isempty(W)
    W=ones(size(V,1),1);
end

if isempty(indStart)
    indStart=1;
end

if ~isfield(optStruct,'numIterationsMax')
    optStruct.numIterationsMax=200;
end

if ~isfield(optStruct,'toleranceLevel')
    optStruct.toleranceLevel=1e-4;
end


if isfield(optStruct,'waitBarOn')
    waitBarOn=optStruct.waitBarOn;
else
    waitBarOn=0;
end


%%

if numSeeds<numel(indStart)
    error('Number of seeds should be equal to or larger than the number of start points')
end

indSeed=indStart;

if waitBarOn
    hw = waitbar(0,'Please wait...');
end

for q=numel(indStart):1:numSeeds
    if waitBarOn
        waitbar(q/numSeeds,hw,['patchMarchDistMapPointSeedIterative ',num2str(round(q/numSeeds*100)),'% complete']);
    end
    
    [D_map] = patchMarchDistMapIterative(F,V,indSeed(1:q),W,optStruct); %New distance map

    if q>=numel(indStart) && q<numSeeds
        [~,indSeed(q+1)]=max(D_map);
    end    
end

if waitBarOn
    close(hw);
end

[~,~,seedIndex]=patchMarchCountPointSeed(F,V,indSeed,numel(indSeed));


% %%
% 
% if numSeeds<numel(indStart)
%     error('Number of seeds should be equal to or larger than the number of start points')
% end
% 
% [E]=patchEdgeLengths(F,V);
% eMin=min(E(:));
% 
% seedIndex=zeros(size(V,1),1);
% indSeed=indStart;
% 
% if waitBarOn
%     hw = waitbar(0,'Please wait...');
% end
% 
% for q=1:1:numSeeds
%     if waitBarOn
%         waitbar(q/numSeeds,hw,['patchMarchDistMapPointSeedIterative ',num2str(round(q/numSeeds*100)),'% complete']);
%     end
%     [D_map] = patchMarchDistMapIterative(F,V,indSeed(1:q),W,optStruct); %New distance map
%     if q==1
%         seedIndex=ones(size(V,1),1)*indSeed(q);
%     else
%         d=(D_map-D_map_now);
%         logicSeedIndex=d<-(eMin/100);
%         seedIndex(logicSeedIndex)=indSeed(q);
%     end
%     if q>=numel(indStart) && q<numSeeds
%         [~,indSeed(q+1)]=max(D_map);
%     end
%     D_map_now=D_map;
% end
% 
% if waitBarOn
%     close(hw);
% end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
