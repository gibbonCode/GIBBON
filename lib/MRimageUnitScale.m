function [Nr,Nr_units]=MRimageUnitScale(N,N_info)

%Scaling data to appropriate units
rescaleSlope=double(N_info.RescaleSlope);
rescaleIntercept=double(N_info.RescaleIntercept);
rescaleType=N_info.RescaleType;

Nr_units=rescaleType;

if isfield(N_info,'ScaleSlope')
    scaleSlope=double(N_info.ScaleSlope);
elseif isfield(N_info,'MRScaleSlope')
    scaleSlope=double(N_info.MRScaleSlope);
else
    scaleSlope=[];
end

if isempty(scaleSlope)
    scaleSlope=1;
    warning('warning, scaleSlope assumed 1!');
end
Nr=((double(N).*rescaleSlope)+rescaleIntercept)./(scaleSlope.*rescaleSlope);



 
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
