function [Mq]=quadrep(M)

% function [Mq]=quadrep(M)
% ------------------------------------------------------------------------
% The input matrix 'M' is mirrored in three directions to create the matrix
% 'Mq'. The size of this matrix is ((2*size(M))-1).
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/09/2008
% ------------------------------------------------------------------------

%% Setting up Mq matrix
Mq=zeros((2*size(M))-1);

%% LEFT SIDE FRONT

%Front top left
Mq(1:size(M,1),1:size(M,2),1:size(M,3))=M;

%Front bottom left
Mq(size(M,1):end,1:size(M,2),1:size(M,3))=flipdim(M,1);

%% RIGHT SIDE FRONT
Mq(1:end,size(M,2):end,1:size(M,3))=flipdim(Mq(1:end,1:size(M,2),1:size(M,3)),2);

%% BACK
Mq(1:end,1:end,size(M,3):end)=flipdim(Mq(1:end,1:end,1:size(M,3)),3);

%% END
 
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
