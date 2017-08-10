function [M]=unpadeven(M_pad,L_even_dim)

num_dims=ndims(M_pad);
siz=size(M_pad)-L_even_dim;
switch num_dims
    case 2
        M=M_pad(1:1:siz(1),1:1:siz(2));
    case 3
        M=M_pad(1:1:siz(1),1:1:siz(2),1:1:siz(3));        
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
