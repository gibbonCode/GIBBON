function c=ivoigtMap(cVoigt)

if ~isvector(cVoigt) %assume that c shall be a 4d array e.g. 4th order tensor
    siz=3*ones(1,4);
elseif isvector(cVoigt) %c shall be a 2D array
    siz=3*ones(1,2); 
    cVoigt(4:end)=cVoigt(4:end)./2;
    secondOrder=1;
else
    %TO DO error catching
end

switch class(cVoigt)
    case 'double'
        c=zeros(siz);
    case 'sym'
        c=sym(zeros(siz));
end

[linearIndexVoigt,linearIndexFull]=tensor2voigtMap(c);
c(linearIndexFull)=cVoigt(linearIndexVoigt);
if secondOrder==1
  [I,J]=ind2sub(siz,linearIndexFull(4:end));
  linearIndexFull2=sub2ind(siz,J,I);
  c(linearIndexFull2)=cVoigt(linearIndexVoigt(4:end));
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
