function D=pathLength(V)

% function [F,V,C]=ind2patch(IND,M,ptype)
% ------------------------------------------------------------------------
%
% This function calculates the "current" curve patch length for each of the
% points of the curve defined by the input argument V. The curve may be
% multidimensional. The output D is a vector of size [size(V,1) 1] whereby
% each entry is defined as the sum of the point-to-point (Euclidean)
% distances (i.e. the curve patch length) leading up to that point. Hence
% it is clear that D(1) is 0 and D(end) is the total curve length. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/05/2013
%------------------------------------------------------------------------

%Compute distance metric
D=zeros(size(V,1),1); %Initialise with zeros (first values stays zero)
D(2:end)=cumsum(sqrt(sum(diff(V,1,1).^2,2)),1);


 
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
