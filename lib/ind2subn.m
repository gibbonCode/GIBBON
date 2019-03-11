function [A] = ind2subn(siz,ind)

% function [A] = ind2subn(siz,ind)
% ------------------------------------------------------------------------ 
% This function is similar to MATLAB's ind2sub function. However the output
% is collected in an array that is numel(ind) x numel(siz). This function
% constructs an array where the columns represent the subscripts indices
% for the input linear indices (ind). 
%
%   See also: ind2sub, sub2ind
% ------------------------------------------------------------------------

%%

siz = double(siz); %The size
ind=double(ind(:)); %The linear indices
numDim=numel(siz); %The number of dimensions
k = cumprod(siz); %Cummulative produc of size

if numDim < 2
    error('Invalid size specified. The number of dimensions should be equal or larger than 2');
end

if any(ind(:)<1) || any(ind>prod(siz))
    error('Index out of range');
end

A=zeros(numel(ind),numDim); %Initializing output array
for q=numDim:-1:1 %For all dimensions
    if q==1 %First 1st dimension
        A(:,1)=rem(ind-1,k(1))+1;        
    else
        p=rem(ind-1,k(q-1)) + 1; %"previous"
        A(:,q)=double((ind - p)/k(q-1) + 1); %Current        
    end    
    ind=p; %Set indices as "previous"      
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
