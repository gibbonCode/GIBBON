function C=makeSymbolicStiffnessTensor(opt)

C=sym(zeros([3 3 3 3]));
I=eye(3,3); 
II1=dyadicProduct(I,I,1);
II3=dyadicProduct(I,I,3);

for i=1:3; 
    for j=1:3;
        for k=1:3; 
            for l=1:3;
                cvar=['c',num2str(i),num2str(j),num2str(k),num2str(l)]; 
         
                
                C(i,j,k,l)=sym(cvar); 
            
            end;
        end; 
    end; 
end;

switch opt
    case 'iso'
        L=(II1+II3)==0;
    case 'transiso'
        
    case'ortho'
        
    case'full'
        L=false(size(C));
    case 'empty'
        L=true(size(C));
end
C(L)=0;

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
