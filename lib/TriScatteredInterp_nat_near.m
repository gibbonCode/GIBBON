function [VI,L]=TriScatteredInterp_nat_near(X,V,XI)

%Natural neighbour interpolation
F=scatteredInterpolant(X,V,'natural');
VI=F(XI);

%Use nearest neighbour for points outside of convex-hull
L=any(isnan(VI),2);
if any(L)
    F=scatteredInterpolant(X,V,'nearest');
    VI(L,:)=F(XI(L,:));
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
