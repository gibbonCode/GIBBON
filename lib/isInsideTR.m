function [varargout]=isInsideTR(varargin)

if nargin<2
    error('Not enough input arguments');
end

TR=varargin{1};
numTR=size(TR.ConnectivityList,1);
QP=varargin{2};

if nargin==3
    ti=varargin{3};
else
    ti=1:numTR;
end
    
%Get barycentric coordinates of points
baryCoords=cartesianToBarycentric(TR,ti(:),QP);

logicFoundEnclosing=all(baryCoords>0,2);

varargout{1}=logicFoundEnclosing;
if nargout==2
    varargout{2}=baryCoords;
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
