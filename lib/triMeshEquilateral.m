function [F,V]=triMeshEquilateral(minV,maxV,pointSpacing)

% function [F,V]=triMeshEquilateral(minV,maxV,pointSpacing)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/28/04
% ------------------------------------------------------------------------
%%

%Region "Field Of View" size
FOV=abs(maxV-minV);

%Calculate number of points in x direction
nx=ceil(FOV(1)./pointSpacing)+1;
xRange=linspace(minV(1)-pointSpacing,maxV(1)+pointSpacing,nx+2);

%Actual point spacing used
pointSpacingReal_X=abs(xRange(1)-xRange(2));

%Calculate number of points in x direction
pointSpacingReal_Y=pointSpacingReal_X.*0.5*sqrt(3);
yRange=minV(2)-pointSpacingReal_Y:pointSpacingReal_Y:maxV(2)+pointSpacingReal_Y;
% yRange=[yRange(1)-pointSpacingReal_Y  yRange  yRange(end)+pointSpacingReal_Y];

%Create normal mesh
[X,Y]=meshgrid(xRange,yRange);

%Offset mesh in X direction to obtain aproximate equilateral triangular mesh
X(1:2:end,:)=X(1:2:end,:)+(pointSpacing/2);

V=[X(:) Y(:)];
V(:,3)=0; 

%Compute connectivity
[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1); %Point indices

%"Forward slash"
TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
F1 = sub2ind(size(X),TRI_I,TRI_J);

%"Backward slash"
TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
F2 = sub2ind(size(X),TRI_I,TRI_J);

F_I=[I(:);I(:)];
logicSlashType=iseven(F_I);

F=[F1(~logicSlashType,:);F2(logicSlashType,:)];

 
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
