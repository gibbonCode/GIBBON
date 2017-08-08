function Cmapped=cmaperise(C,cmap,clim)

C=C(:);
p=(C-clim(1))./(clim(2)-clim(1));
p(p<0)=0;
p(p>1)=1;
IND_cmap=round((p*(size(cmap,1)-1))+1);
Cmapped=cmap(IND_cmap,:);

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
