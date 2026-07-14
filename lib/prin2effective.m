function [D_eff]=prin2effective(varargin)

% function [D_eff]=prin2effective(D1,D2,D3,typeFlag)
% ------------------------------------------------------------------------
% This function computes the effective/Von Mises stress or strain based on
% the input principal components. If typeFlag='stress' then the von Mises
% stress is computer. If typeflag='strain' the effective strain is
% computed. 
% 
% 2023/05/11: Created
% 2023/05/30: Added documentation/description
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        D=varargin{1};
        typeFlag=varargin{2};
        D1=D(:,1,:);
        D2=D(:,2,:);
        D3=D(:,3,:);        
    case 4
        D1=varargin{1};
        D2=varargin{2};
        D3=varargin{3};  
        typeFlag=varargin{4};
    otherwise
        error('Wrong number of input arguments.');
end

%%

%Compute Von Mises / effective metric
switch typeFlag
    case 'stress' %Von Mises Stress
        D_eff =         sqrt( ((D1-D2).^2 + (D2-D3).^2 + (D1-D3).^2)./2 );
    case 'strain' %Effective strain
        D_eff = (1/3).* sqrt( ((D1-D2).^2 + (D2-D3).^2 + (D1-D3).^2).*2 );
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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
