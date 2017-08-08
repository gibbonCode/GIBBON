function [E,V]=rhombicDodecahedronHex(r)

%% Get rhombic dodecahedron

[F,V]=rhombicDodecahedron(r);

%% Add centre point

V(end+1,:)=zeros(1,3);

%% Define the elements

E1=[fliplr(F(1,:)) fliplr([size(V,1) F(5,3) F(5,4) F(8,1)]) ];
E2=[fliplr(F(2,:))  fliplr([F(9,2) F(10,3) size(V,1) F(9,1)]) ];
E3=[fliplr(F(3,:)) fliplr([size(V,1)  F(7,3) F(6,4) F(6,1)]) ];
E4=[fliplr((F(4,:))) fliplr([F(12,2) F(12,3) size(V,1) F(11,1)])  ];
E=[E1;E2;E3;E4];

 
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
