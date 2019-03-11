function [X]=ppInv(varargin)

% function [X]=ppInv(pp,Y,outputOpt)
% ------------------------------------------------------------------------
%
%
%
% Based on PPInverse from:
% http://web.mit.edu/~paul_s/www/14.170/matlab/dp/v5-soln/approx/ppcreate.m
% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-03-23, 2007-03-01
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/08/31 Updated for GIBBON, expanded for vector input
%------------------------------------------------------------------------

%%

pp=varargin{1};
Y=varargin{2};

switch nargin
    case 2
        outputOpt='vec';
    case 3
        outputOpt=varargin{3};
end

dxbr=diff(pp.breaks);
tol=100*eps;

switch outputOpt
    case 'cell'
        X=cell(size(Y));
    case 'vec'
        X=nan(size(Y));
end

for q=1:1:numel(Y)
    x=[];
    y=Y(q);
    ppn=pp;
    ppn.coefs(:,end)=pp.coefs(:,end)-y;  % shift all polys by y
    
    for k=1:ppn.pieces
        p=ppn.coefs(k,:);     % k-th poly to investigate
        inz=find(p);         % index of nonzero elements of poly
        fnz=inz(1);          % first nonzero element
        lnz=inz(end);        % last nonzero element
        
        p=p(fnz+1:lnz)/p(fnz);        % strip leading,trailing zeros, make monic
        r=zeros(1,ppn.order-lnz);      % roots at zero
        
        if (lnz-fnz)>1                   % add nonzero roots
            a=diag(ones(1,lnz-fnz-1),-1); % form companion matrix
            a(1,:)=-p;
            r = [r eig(a).'];             % find eigenvalues to get roots
        end
        
        i=find(abs(imag(r))<tol & ...          % find real roots with
            real(r)>=0 & ...                % nonnegative real parts that are
            real(r)<dxbr(k)-eps(dxbr(k)));  % less than the next breakpoint
        if ~isempty(i)
            x=[x; ppn.breaks(k)+r(i)];
        end
    end
    
    if ~isempty(x)
        switch outputOpt
            case 'cell'
                X{q}=x;
            case 'vec'
                if numel(x)>1
                    error('Multiple solutions exist. Try changing output type to: cell')
                else
                    X(q)=x;
                end
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
