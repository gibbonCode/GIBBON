function [varargout]=patchPointDist(V,F,Vp,dL)

methodOpt='near-norm';

switch methodOpt
    case 'near-norm'
        
        [Nn,Vn]=patchNormal(F,V); %Get face normals and normal vector origins
        [Dn,indMin]=minDist(Vp,Vn); %Get distances to origins and find closest faces
        
        %Attempt to find orthogonal projection to closest face
        N=Nn(indMin,:); %Order the normals
        W=Vp-Vn(indMin,:);
        Dm=(dot(N,W,2));
        D=abs(Dm);        
        Dp=Dm(:,ones(1,3)).*N;
        Vpd=Vp-Dp;
        
        %Test if projection is valid
        TR = triangulation(F,V);
        B = cartesianToBarycentric(TR,indMin,Vpd);        
        logicOutside=any(B<-dL,2);
        
        %Replace invalid projections with face centre points
        D(logicOutside)=Dn(logicOutside); %Replace distances        
        Dp(logicOutside,:)=Vp(logicOutside,:)-Vn(indMin(logicOutside),:); %Replace difference vectors         
        
    otherwise
        error('Wrong method option chosen');
end

switch nargout
    case 1
        varargout{1}=D;
    case 2
        varargout{1}=D;
        varargout{2}=Dp;    
    case 3
        varargout{1}=D;
        varargout{2}=Dp;
        varargout{3}=logicOutside;
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
