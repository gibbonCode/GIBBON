function [viewProfile]=getViewProfile

% function [viewProfile]=getViewProfile
% ------------------------------------------------------------------------
% This function checks for the current view profile setting. If no setting
% is made then the default 'CAD' option is returned.
%
% See also: setViewProfile, vcw
% 
% 2023/05/30 Created by Kevin Moerman
% ------------------------------------------------------------------------

%%

settingSet=settings; %Get settings
viewProfile='CAD'; %Initialize as default CAD
if hasGroup(settingSet,'GIBBON') %Prior GIBBON settings group exists -> Check for setting
    if hasSetting(settingSet.GIBBON,'vcw_profile') %Prior vcw_profile setting exists -> Update setting
        viewProfile=settingSet.GIBBON.vcw_profile.PersonalValue; %Access set value
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