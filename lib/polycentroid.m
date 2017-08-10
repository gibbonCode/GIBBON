function [Xc,Yc]=polycentroid(X,Y)

%N.B. 
% Assumes row vectors or matrices whereby each row describes a polygon with
% points appearing in the order defining the polygon

meanX=mean(X,2);
meanY=mean(Y,2);
X=X-meanX*ones(1,size(X,2));
Y=Y-meanY*ones(1,size(Y,2));

A = polyarea(X',Y');
Xc=sum((X(:,1:end-1)+X(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));
Yc=sum((Y(:,1:end-1)+Y(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));

Xc=Xc+meanX;
Yc=Yc+meanY;

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
