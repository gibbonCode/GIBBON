function [varargout]=minPolyTwist(V1,V2)

%Get number of points for curves
n1=size(V1,1);
n2=size(V2,1);

if n1==n2
    
    %Fix overall direction if required
    [L1]=isPolyClockwise(V1);
    [L2]=isPolyClockwise(V2);
    flipReq=L1~=L2;
    
    if flipReq
        flipFlag=1;
        V2t=flipud(V2);
    else
        flipFlag=0;
        V2t=V2;
    end
    
    %Find point order for minimum torsion
    SSQD=zeros(n1,1);
    for q=1:1:n1
        if q>1
            indSort=[q:n1 1:q-1];
        else
            indSort=1:n1;
        end
        V2s=V2t(indSort,:);
        SSQD(q)=sum(sqrt(sum((V1-V2s).^2,2)).^2);
    end
    
    %Get sort order
    [~,qMin]=min(SSQD);
    if qMin>1
        indSort=[qMin:n1 1:qMin-1];
    else
        indSort=1:n1;
    end
    
    %Fix order if flipping was required
    if flipReq
        indSort=n1-(indSort-1);
    end
    
    %Reorder points
    V2f=V2(indSort,:);
    
    switch nargout
        case 1
            varargout{1}=V2f;
        case 2
            varargout{1}=V2f;
            varargout{2}=indSort;
        case 3
            varargout{1}=V2f;
            varargout{2}=indSort;
            varargout{3}=flipFlag;            
    end
else
    error('size(V1,1) should equal size(V2,1)');
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
