function M=checkerBoard3D(varargin)

% function M=checkerBoard3D(siz)
% ------------------------------------------------------------------------
% This function creates a checkboard image of the size siz whereby elements
% are either black (0) or white (1). The first element is white.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% 2018/12/10 Added checkboard block size input
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        siz=varargin{1};
        blockSize=1;
    case 2        
        siz=varargin{1};        
        blockSize=varargin{2};
end

if ~isrounded(blockSize) || blockSize<1
    error('Block size should be positive integer');
end

%Coping with 1D or 2D input
if numel(siz)==2
    siz(3)=1; 
elseif numel(siz)==1
    siz(2:3)=[1 1];
end

%%

if blockSize==1
    [I,J,K]=ndgrid(1:1:siz(1),1:1:siz(2),1:1:siz(3));
else    
    [I,J,K]=ndgrid(indexRange(blockSize,siz(1)),indexRange(blockSize,siz(2)),indexRange(blockSize,siz(3)));
end
logic_ij=((iseven(I)| iseven(J)) & ((iseven(I)~=iseven(J))));
M=false(siz);
M(iseven(K))=logic_ij(iseven(K));
M(~iseven(K))=~logic_ij(~iseven(K));

end

%%
function i=indexRange(blockSize,s)

i=repmat(1:blockSize:s,blockSize,1); 
i=i(:);
if numel(i)<s
    ii=i;
    i=siz(1)*ones(1,s);
    i(1:num_i)=ii;
else
    i=i(1:s);
end
[~,~,i]=unique(i);
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
