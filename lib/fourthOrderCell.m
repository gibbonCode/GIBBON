function [CM]=fourthOrderCell(C)

CM=cell(3,3);

for i=1:1:3
    for j=1:1:3
        switch class(C)
            case 'sym'
                C_sub=sym(zeros(3,3));
            otherwise
                C_sub=zeros(3,3);                
        end
        
        for k=1:1:3
            for l=1:1:3                
                C_sub(k,l)=C(i,j,k,l);
            end
        end
        CM{i,j}=C_sub;
    end
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
