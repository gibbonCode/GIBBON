function [xb]=boxconstrain(varargin)

% function [xb]=boxconstrain(x,lb,ub,m,s)
%-------------------------------------------------------------------------
% The boxconstrain function can be used to constrain parameters from
% [-inf,inf] to the range [lb ub]. The sigmoidal tanh function is used to
% do this mapping. The middle of the sigmoid, denoted by the parameter m
% need not be (lb+ub)/2. If x=m xb=m. 
% The optional parameter s sets the slope of the sigmoid at
% x=m, the default slope is 1. 
% 
% 2018/06/22 Created to replace parLimNat
%-------------------------------------------------------------------------

%%
% Parse input

switch nargin    
    case 3
        x=varargin{1};
        lb=varargin{2};
        ub=varargin{3};
        m=(ub+lb)/2;
        s=1;
    case 4
        x=varargin{1};
        lb=varargin{2};
        ub=varargin{3};
        m=varargin{4};
        s=1;
    case 5
        x=varargin{1};
        lb=varargin{2};
        ub=varargin{3};
        m=varargin{4};
        s=varargin{5};        
end

%%
% Check for bad parameters
if s<=0
    error('The slope should larger than 0');
end

if lb>ub
    error('Lower limit exceeds upper limit');
end

if ub<lb
    error('Upper limit smaller than lower limit');
end

if abs(ub-lb)<eps
    error('Limits are too similar');
end

if m<lb || m>ub
   error('Middle should be between limits'); 
end

if m==lb || m==ub
   warning('Middle coincides with a boundary, transition will not be smooth'); 
end

%%
% Compute mapping to constrained space based on tanh function

if isempty(x)
    xb=[];
else
    xb=(heaviside(m-x).*(m+tanh((x-m)./(m-lb).*s).*(m-lb)))... %part on left of middle
      +(heaviside(x-m).*(m+tanh((x-m)./(ub-m).*s).*(ub-m)));   %part on right of middle
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
