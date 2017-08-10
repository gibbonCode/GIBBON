function [C]=linspacen(A,B,n)

% function [C]=linspacen(A,B,n)
% ------------------------------------------------------------------------
% This function is a generalization of the linspace function to N
% dimensions. The output C is a matrix of size [size(A) n] such that "it
% goes from A to B in n steps in the last dimention. The input variables A
% and B (scalars, vectors or matrices). For scalar input this function is
% equivalent to linspace (but slower due to repmat operation).
% Clearly the inputs A and B should have the same size.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% Updated: 15/07/2010
%------------------------------------------------------------------------

%%

v=isvector(A);
size_A=size(A);
A=A(:);
B=B(:);

C=repmat(A,[1,n])+((B-A)./(n-1))*(0:1:n-1);
% C(:,end)=B;

if v~=1
    C=reshape(C,[size_A n]);
end


% %Alternative without REPMAT
% q=permute(linspace(0,1,n),(numel(size(A))+2):-1:1);
% Q=zeros([size(A) n]);
% Q=q(ones(size(A,1),1),ones(size(A,2),1),ones(size(A,3),1),:);
% A=A(:,:,:,ones(n,1));
% B=B(:,:,:,ones(n,1));
% D=B-A;
% D=D(:,:,:,ones(n,1));
% C=A+Q.*D;


 
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
