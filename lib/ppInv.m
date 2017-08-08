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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
