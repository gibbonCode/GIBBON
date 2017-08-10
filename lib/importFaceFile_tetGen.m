function [varargout]=importFaceFile_tetGen(fileName)

fid=fopen(fileName,'r');
[A]=textscan(fid,'%d %d %d %d %f','HeaderLines',1,'Delimiter',' ','CommentStyle','Shell','MultipleDelimsAsOne',1);
fclose(fid);
faceID=double(A{1});
F=nan(max(faceID),3);
F(faceID,:)=double([A{2} A{3} A{4}]);
faceBoundaryID=double(A{5});

if all(isnan(faceBoundaryID(:)))
    faceBoundaryID=-ones(size(F,1),1);
end

varargout{1}=faceID;
varargout{2}=F;
varargout{3}=faceBoundaryID;

 
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
