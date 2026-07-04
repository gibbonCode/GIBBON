function cMap=spectral(varargin)

%%

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%

cMap=[0.368627450980392,0.309803921568627,0.635294117647059;...
      0.196078431372549,0.533333333333333,0.741176470588235;...
      0.400000000000000,0.760784313725490,0.647058823529412;...
      0.670588235294118,0.866666666666667,0.643137254901961;...
      0.901960784313726,0.960784313725490,0.596078431372549;...
      1,1,0.749019607843137;...
      0.996078431372549,0.878431372549020,0.545098039215686;...
      0.992156862745098,0.682352941176471,0.380392156862745;...
      0.956862745098039,0.427450980392157,0.262745098039216;...
      0.835294117647059,0.243137254901961,0.309803921568627;...
      0.619607843137255,0.00392156862745098,0.258823529411765];

if n~=size(cMap,1)
    cMap=resampleColormap(cMap,n);
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
