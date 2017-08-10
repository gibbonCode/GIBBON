function [AB]=vectorTensorProductArray(a,b)

%% 

%Determine multiplication order
if size(a,2)==3 %type B*A
   B=a;
   A=b;
   multType='post';
else %type A*B
    B=b;
    A=a; 
    multType='pre';
end

%Perform product
switch multType
    case 'pre' %type A*B
        %  A1_1*B1 + A1_2*B2 + A1_3*B3
        %  A2_1*B1 + A2_2*B2 + A2_3*B3
        %  A3_1*B1 + A3_2*B2 + A3_3*B3
        AB=[A(:,1).*B(:,1)+A(:,4).*B(:,2)+A(:,7).*B(:,3) ...
            A(:,2).*B(:,1)+A(:,5).*B(:,2)+A(:,8).*B(:,3) ...
            A(:,3).*B(:,1)+A(:,6).*B(:,2)+A(:,9).*B(:,3)];
    case 'post' %type B*A
        %  A1_1*B1 + A1_2*B2 + A1_3*B3
        %  A2_1*B1 + A2_2*B2 + A2_3*B3
        %  A3_1*B1 + A3_2*B2 + A3_3*B3
        AB=[A(:,1).*B(:,1)+A(:,2).*B(:,2)+A(:,3).*B(:,3) ...
            A(:,4).*B(:,1)+A(:,5).*B(:,2)+A(:,6).*B(:,3) ...
            A(:,7).*B(:,1)+A(:,8).*B(:,2)+A(:,9).*B(:,3)];
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
