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
