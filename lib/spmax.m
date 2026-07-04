function [varargout]=spmax(varargin)

% function [maxVal,maxInd]=spmax(A,B,vecdim,nanflag,logicRelevant,nanOut)
% ------------------------------------------------------------------------
%
%
% Kevin Moerman
% ------------------------------------------------------------------------

%%

A=varargin{1};

if nargin>=3
    dimDir=varargin{3};
else
    dimDir=1;
end

switch nargin
    case 5        
        logicRelevant=varargin{5};
        varargin=varargin(1:4);
        nanOut=0;
        if isempty(varargin{4})
            varargin=varargin(1:3);
        end
    case 6
        logicRelevant=varargin{5};
        nanOut=varargin{6};
        varargin=varargin(1:4);
        if isempty(varargin{4})
            varargin=varargin(1:3);
        end
    otherwise
        logicRelevant=[];
        nanOut=0;
end

if isempty(logicRelevant)
    logicRelevant=A~=0;
end

logicRelevantRow=full(any(logicRelevant,dimDir));

%%

minOffset=min(A(logicRelevant));

A(logicRelevant)=A(logicRelevant)-minOffset; %Shift to negative to zeros loose
varargin{1}=A;
[maxVal,maxInd]=max(varargin{:});
maxVal(logicRelevantRow)=maxVal(logicRelevantRow)+minOffset; %Shift back
if nanOut==1
    maxVal(~logicRelevantRow)=NaN;
end

%% Collect output
varargout{1}=maxVal;
varargout{2}=maxInd;

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
