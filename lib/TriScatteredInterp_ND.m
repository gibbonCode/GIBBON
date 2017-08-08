function [DI]=TriScatteredInterp_ND(DT,D,XI,InterpMethod)

if ~isa(DT,'DelaunayTri') %if DT is not a delaunay tesselation
    DT=delaunayTriangulation(DT); %assuming DT are coordinates, replace by Delaunay tesselation
end

DI=nan(size(XI,1),size(D,2)); %Allocate DI
for q=1:size(D,2)% loop over dimensions
    switch InterpMethod
        case 'nat_near' %natural in chull, neirest outside chull
            [DI(:,q),~]=TriScatteredInterp_nat_near(DT,D(:,q),XI);
        otherwise %TriScatterInterp can handle other methods
            F = scatteredInterpolant(DT,D(:,q),InterpMethod);
            %Construct interpolator
            DI(:,q)=F(XI); %Interpolate
    end
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
