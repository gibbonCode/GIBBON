function axLim=axisLim(varargin)

% function axLim=axisLim(V)
% ------------------------------------------------------------------------
% Computes tight axis limits for the input coordinates V. The coordinates
% may be a 3-dimensional array where the 3rd direction could reflect
% coordinates as a function of time for instance. 
% The function also supports multiple inputs e.g.:
%
% axLim=axisLim(V1,V2,V3)
%
% In this case the limits are based on all sets simultaneously. Note some
% may be 3-dimensional arrays (e.g. nx3xm) while others may be 2-dimensional arrays. 
%
% 2022/11/08 Updated to allow for multiple coordinate set input
% ------------------------------------------------------------------------
%%

for q=1:1:nargin %Loop over all coordinate sets

    V=varargin{q};

    try
        minV=min(V,[],[1 3]);
        maxV=max(V,[],[1 3]);
    catch
        minV=min(min(V,[],1),[],3);
        maxV=max(max(V,[],1),[],3); 
    end    
    
    if q>1 %Compare to previous                
        %Check for dimensionality compatibility
        if numel(minVn)==2 && numel(minV)==3
            minVn(3)=0;
            maxVn(3)=0;
        elseif numel(minV)==2 && numel(minVn)==3
            minV(3)=0;
            maxV(3)=0;
        end
        minV=min(minV,minVn);
        maxV=max(maxV,maxVn);
    end
    
    %Keep track of previous
    minVn=minV;
    maxVn=maxV;
end

logicSame= abs(minV-maxV) < eps(abs(minV-maxV));
minV(logicSame)=minV(logicSame)-1;
maxV(logicSame)=maxV(logicSame)+1;

%Create axis limit format
axLim=[minV(:) maxV(:)]';
axLim=axLim(:)'; %Format is X X Y Y Z Z

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
