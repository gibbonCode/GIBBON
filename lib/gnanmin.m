function [varargout]=gnanmin(varargin)
% function [y,ind]=gnanmin(A,[],dimDir)
%------------------------------------------------------------------------
% this function is a replacement for nanmin which is part of the
% statistics and machine learning toolbox. 
% The use of min with the omitnan flag is used or if this does not work
% (for old MATLAB versions) a custom implementation is used. 
%
% Change log:
% 2020/01/09 Created
%------------------------------------------------------------------------

%%
if nargin==1
    varargin{2}=[];
end

try 
    %Use 'omitnan' flag    
    if nargin==1
        varargin{2}=[];
    end
    if nargout==2
        [y,ind] = min(varargin{:},'omitnan');
    else
        y = min(varargin{:},'omitnan');
    end
catch
    %Use custom version     
    A=varargin{1}; %First array
    B=varargin{2}; %Second array or number or empty
    
    if nargin==3
        dimDir=varargin{3};
    else
        dimDir=1;
    end
    
    LA=isnan(A);    
    maxVal=max(A(~LA));
    
    if ~isempty(B)
        LB=isnan(B);
        maxVal=maxVal+max(B(~LB));
    end
    if isempty(maxVal)
        maxVal=0;
    end
    
    AV=A;
    AV(LA)=1+maxVal;    
    if ~isempty(B)
        if nargout==2
            error('If two matrices are compared two output arguments are not supported.');
        end
        BV=B;
        BV(LB)=1+maxVal;                
        y=min(AV,BV); 
        y(LA&LB)=NaN;
        y(~LA&LB)=A(~LA&LB);
        y(LA&~LB)=B(LA&~LB);        
    else
        [y,ind]=min(AV,[],dimDir);            
        y(all(LA,dimDir))=NaN; %replace all NaN directions
    end
    
end

varargout{1}=y;
if nargout==2
    varargout{2}=ind;
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
