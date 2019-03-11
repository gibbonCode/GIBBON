function [B]=gcombvec(varargin)

numElements=cellfun(@numel,varargin);
numVecs=numel(varargin);
siz=prod(numElements);
B=zeros(numVecs,siz);

for q=1:1:numVecs
     vecNow=varargin{q};
     if q==1
         nRep1=1;
     else
         nRep1=prod(numElements(1:q-1));
     end     
     if q==numVecs
         nRep2=1;
     else
         nRep2=prod(numElements(q+1:end));         
     end
     w=repmat(vecNow,[nRep1 1]);    
     w=repmat(w(:),[nRep2 1]);          
     B(q,:)=w;
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
