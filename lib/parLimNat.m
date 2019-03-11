function [varargout]=parLimNat(varargin)

% function [xx,S]=parLimNat(xx_c,[xx_min xx_max],x);
% ------------------------------------------------------------------------
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
%
% 2015/06/30 Updated for GIBBON, fixed case for limits equal to centre
% 2015/08/31 Fixed behaviour for empty x, changed to varargin
% 2018/06/21 Increased sampling to make output smoother
% 2018/06/22 Added warning to refer to boxconstrain
%------------------------------------------------------------------------

%%

warning('This function is obsolete and will be removed in future releases, use the boxconstrain function instead');

%%

xx_c=varargin{1};
xxlim=varargin{2};

switch nargin       
    case 2
        x=[];
    case 3 
        x=varargin{3};
end

%%

if isempty(xx_c)
    xx_c=mean(xxlim);
end

if diff(xxlim)<eps(xx_c)
    error('Limits xxlim are too similar');
end

if xxlim(1)>xx_c
    error('Lower limit should be lower or equal to centre');
end

if xxlim(2)<xx_c
    error('Upper limit should be higher or equal to centre');
end

%%

wn=abs(xx_c-xxlim(1));
L_wn=wn<eps(xx_c);
if L_wn
   warning('Centre should not coincide with lower limit. Transition will not be smooth'); 
end

wp=abs(xxlim(2)-xx_c);
L_wp=wp<eps(xx_c);
if L_wp
   warning('Centre should not coincide with upper limit. Transition will not be smooth'); 
end

%%

nc=100; 

if L_wn
    xtn=xx_c;
else
    xtn=linspace(xx_c-2*wn,xx_c,nc);
end

if L_wp
    xtp=xx_c;
else
    xtp=linspace(xx_c,xx_c+2*wp,nc);
end

xt=[xtn xtp(2:end)];
xtlim=[min(xt(:)) max(xt(:))];

xxt=xt;

Ln=xt<=xx_c;
if L_wn
    xxt(Ln)=0;
else
    xxt(Ln)=erf(1/wn*2*(xt(Ln)-xx_c));
    xxt(Ln)=xxt(Ln)*wn;
end

Lp=xt>xx_c;
if  L_wp
    xxt(Lp)=[];
else
    xxt(Lp)=erf(1/wp*2*(xt(Lp)-xx_c));
    xxt(Lp)=xxt(Lp)*wp;
end
xxt=xxt+xx_c;

%% Construct and evaluate Piecewise Cubic Hermite Interpolating Polynomial 
pp = pchip(xt,xxt); 
if ~isempty(x)
    xx=ppval_extrapVal(pp,x,xtlim,xxlim);
else
    xx=[];
end

%% Create output
varargout{1}=xx; 
switch nargout
    case 2      
        %Create additional output structure
        S.pp=pp;
        S.interplim=xtlim;
        S.extrapVal=xxlim;        
       
        varargout{2}=S;
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
