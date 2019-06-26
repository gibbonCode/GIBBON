function [varargout]=axis_inc(varargin)

% function [h]=axis_inc(s,ha)
% ------------------------------------------------------------------------
% This function widens the axis limits for the axis with the handle |hs|
% using the scale factor |s| (e.g. if s=2 the axis width is doubled). The
% scaling parameter |s| can be a single scalar or a vector such that each
% axis direction is scaled differently.  
% 
% Change log: 
% 2019/06/26 Added variable input/output handling
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        s=varargin{1};
        ha=gca; 
    case 2
        s=varargin{1};
        ha=varargin{2};        
end

if numel(s)==1
    s=s*ones(1,3);    
end

%% Set new axis limits

%Get axis limits
limSet={'XLim','YLim','ZLim'};

for q=1:1:numel(limSet)
    try
        limNow=get(ha,limSet{q});
        w=limNow(2)-limNow(1);
        if w<=eps(0)
            w=1; %Force 1 if too small
        end
        limInc=w*(s(q)-1)/2; %Axis limit increment
        limNew=limNow+[-limInc limInc];
        set(ha,limSet{q},limNew);
    catch
        
    end
end

%% Collect output
varargout{1}=ha;

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
