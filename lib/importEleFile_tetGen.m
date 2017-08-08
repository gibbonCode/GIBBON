function [varargout]=importEleFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
elementID=double(A{1});
E=nan(max(elementID),4);
E(elementID,:)=double([A{2} A{3} A{4} A{5}]);
elementMaterialID=double(A{6});

if all(isnan(elementMaterialID(:)))
    elementMaterialID=-ones(size(E,1),1);
end

varargout{1}=elementID;
varargout{2}=E;
varargout{3}=elementMaterialID;

 
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
