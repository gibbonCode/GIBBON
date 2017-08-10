function [varargout]=euler2DCM(E)

% ----------------------------------------------------------------------
% function [R,Ri]=euler2DCM(E)
% 
% This function generates a rotation matrices based on the Euler angles in
% E.
%
% Kevin Mattheus Moerman, 01/06/2011
% kevinmoerman@hotmail.com
% ----------------------------------------------------------------------


switch class(E)
    case 'double'
        R=zeros(3,3,size(E,1));
    case 'sym'
        R=sym(zeros(3,3,size(E,1)));
end
Ri=R;

for q=1:1:size(E,1)
    
    Rx=[1        0              0;...
        0        cos(E(q,1))  -sin(E(q,1));...
        0        sin(E(q,1))   cos(E(q,1))];
    
    Ry=[cos(E(q,2))  0        sin(E(q,2));...
            0        1        0;...
        -sin(E(q,2)) 0        cos(E(q,2))];
    
    Rz=[cos(E(q,3))  -sin(E(q,3)) 0;...
        sin(E(q,3))  cos(E(q,3))  0;...
        0        0        1];
    
    Rxyz=Rx*Ry*Rz;
    R(:,:,q)=Rxyz;
    
    Ri(:,:,q)=inv(Rxyz);
end

varargout{1}=R; 
varargout{2}=Ri; 

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
