function [Ri]=interp_polar(TZ,radiusData,TZi,interp_method)

if size(TZ,2)==1 %interp1 function
    thetaData=TZ;
    thetaData=[(thetaData-2.*pi); thetaData; (thetaData+2.*pi)];
    
    thetaData=thetaData.*min(radiusData(:));
    TZi=TZi.*min(radiusData(:));
    
    radiusData=repmat(radiusData,[3,1]);
    if ~ischar(interp_method);
        Ri = csaps(thetaData,radiusData,interp_method,TZi);
    else
        Ri=interp1(thetaData,radiusData,TZi,interp_method);
    end
elseif size(TZ,2)==2 
    thetaData=[(TZ(:,1)-2.*pi); TZ(:,1); (TZ(:,1)+2.*pi)];
    radiusData=repmat(radiusData,[3,1]);
    ZData=repmat(TZ(:,2),[3,1]);
    
    %Scaling angular units to length units using min(R(:))
    thetaData=thetaData.*min(radiusData(:));
    TZi(:,1)=TZi(:,1).*min(radiusData(:));
    
    %Interpolating
    if strcmp(interp_method,'natural') || strcmp(interp_method,'linear') || strcmp(interp_method,'nearest') %TriScatterdInterp function
        F=scatteredInterpolant([thetaData ZData],radiusData,interp_method);
        Ri=F(TZi);
    else %Griddata function
        Ri = griddata(thetaData,ZData,radiusData,TZi(:,1),TZi(:,2),interp_method);
    end
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
