function setViewProfile(profileName)

% function setViewProfile(profileName)
% -----------------------------------------------------------------------
% This function sets the view profile used by vcw. 
% 
% See also: getViewProfile, vcw
% 
% 2023/05/25 Created by Kevin Moerman
% ------------------------------------------------------------------------

%% Parse input

validSet={'default','CAD','febio','touchpad'};
if ~ismember(profileName,validSet)    
    error('Invalid profile type provided'); 
end

%% Create vcw setting
s=settings;

if hasGroup(s,'GIBBON') %Prior GIBBON settings group exists -> Check for setting
    if ~hasSetting(s.GIBBON,'vcw_profile') %Prior vcw_profile setting exists -> Update setting        
        addSetting(s.GIBBON,'vcw_profile'); % Add vcw_profile setting
    end
else %No GIBBON settings at all -> Make group and setting
    addGroup(s,'GIBBON'); % Add GIBBON settings group
    addSetting(s.GIBBON,'vcw_profile'); % Add vcw_profile setting    
end
s.GIBBON.vcw_profile.PersonalValue=profileName; %Make setting as personal setting

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
