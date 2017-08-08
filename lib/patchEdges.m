function varargout=patchEdges(varargin)

% function [E,indUni1,indUni2]=patchEdges(F,uniOpt)
% -----------------------------------------------------------------------
% E=patchEdges(F,uniOpt)
% Uses the input faces array F to compute an edge array E. If uniOpt==1
% then the output array contain unique edges (irrespective of node order
% such that e.g. [4 1] and [1 4] are seen as the same edge). If uniOpt~=1
% then the double edges are maintained. If only one input is provided it is
% assumed to represent F and the default, uniOpt=0, is used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/17
% 2016/01/19 Added additional output in relation to unique operation
%------------------------------------------------------------------------

%% PARSE INPUT
F=varargin{1};
switch nargin
    case 1
        uniOpt=0;
    case 2
        uniOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end
%% DERIVE NON-UNIQUE EDGES MATRIX
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

%% REMOVE DOUBLE ENTRIES IF DESIRED

if uniOpt==1
    Es=sort(E,2); %Sorted so [1 4] and [4 1] are seen as the same edge
    [~,indUni1,indUni2]=unique(Es,'rows'); %Get indices for unique edges
    E=E(indUni1,:);    
end

%% Collect output

varargout{1}=E; 

if nargout>1 && uniOpt==1
    varargout{2}=indUni1;
    varargout{3}=indUni2;
elseif nargout>1 && uniOpt~=1
    error('Multiple outputs only available if uniOpt=1');
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
