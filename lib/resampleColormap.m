function [cMap_i]=resampleColormap(cMap,n)

% function [cMap_i]=resampleColormap(cMap,n)
% ------------------------------------------------------------------------
%
% This function resamples the input colormap cmap using n steps. Resampling
% is based on linear interpolation.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log:
% 2014/09/25
% 2018/06/13 KMM: Added handling of single color colormap
% 2023/906/06 KMM: Added L,a,b colorspace based interpolation 
%------------------------------------------------------------------------

%%

%Get range
maxC=max(cMap(:));
minC=min(cMap(:));

%Resample colormap
if n~=size(cMap,1) %If resampling is required

    % Attempt conversion to L,a,b color space for homogeneous interpolation
    try
        cMap=rgb2lab(cMap); %Convert RGB to L,a,b color space
    catch
        %waring('rgb2lab and/or lab2rgb not found, using RGB based interpolation ')
    end

    if size(cMap,1)==1        
        cMap_i=cMap(ones(n,1),:);
    else
        ind=(1:1:size(cMap,1))';
        ind_i=linspace(1,size(cMap,1),n)';
        cMap_i=zeros(n,size(cMap,2));
        
        %Interpolate color data
        for q=1:1:size(cMap,2)
            cMap_i(:,q)=interp1(ind,cMap(:,q),ind_i,'linear');
        end
    end

    %Attempt conversion from L,a,b to RGB
    try
        cMap_i=lab2rgb(cMap_i); %Convert back to RGB from L,a,b color space

        %Fix minor overshoot (e.g. due to numerical precission)
        cMap_i(cMap_i>maxC)=maxC;
        cMap_i(cMap_i<minC)=minC;
    catch
        %waring('rgb2lab and/or lab2rgb not found, using RGB based interpolation ')
    end
else
    cMap_i=cMap;
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
