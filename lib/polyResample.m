function [Vi]=polyResample(V,d_n,interpOpt,interpMethod)

%Parameterize curve using curve length
D = polyCurveLength(V); %Create distance space

if nnz(D>max(eps(V),[],2))<(size(D,1)-1)
    warning('Resampling not possible due to non-unique points (zero distances detected)');
    Vi=[];
else
    
    %Creating even space steps allong curve
    switch interpOpt
        case 1 %d_n is taken to be number of desired points
            D_geo=linspace(0,D(end),d_n);
        case 2 %d_n is taken to be desired point spacing
            d_n_fix=D(end)/round(D(end)/d_n);
            if d_n_fix~=d_n
                disp(['Point spacing was set to: ', num2str(d_n_fix)])
            end
            D_geo=0:d_n_fix:D(end);
    end
    
    %Interpolating curve points
    Vi=zeros(numel(D_geo),size(V,2)); %Initialising Vi
    for qd=1:1:size(V,2)
        Vi(:,qd) = interp1(D,V(:,qd),D_geo,interpMethod);
    end
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
