function [Xs,Ys,Zs]=snap2grid(X,Y,Z,c,v)

Xs=c(1)+round((X-c(1))./v(1)).*v(1);
Ys=c(2)+round((Y-c(2))./v(2)).*v(2);
Zs=c(3)+round((Z-c(3))./v(3)).*v(3);

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
