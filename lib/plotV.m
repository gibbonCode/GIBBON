function varargout=plotV(V,varargin)

% function varargout=plotV(V,varargin)
% ------------------------------------------------------------------------
% This function is just like plot3 except instead of having to specify the
% X, Y, and Z coordinates the user can provide a single Nx3 array V for the
% coordinates. The optional input varargin may contain additional inputs as
% possible for plot3. 
% 
% See also: plot3
%
% Change log: 
% 2021/07/22 KMM Change empty output to be empty graphics object rather
% than []. Switched to more compact/efficient code. Added comments. 
% 
% ------------------------------------------------------------------------

%%

if isa(V,'cell') %Recursively loop over cell entries
    h=gobjects(numel(V));
    for q=1:1:numel(V)                
        h(q)=plotV(V{q},varargin{:});
    end
    varargout{1}=h;
else
    if ~isempty(V)
        if size(V,2)==2 %Force 3D
            V(:,3)=0; %Add zeros for Z if input is 2D
        end
        if nargout==1
            varargout{1}=plot3(V(:,1),V(:,2),V(:,3),varargin{:}); %Add handle to output
        else
            plot3(V(:,1),V(:,2),V(:,3),varargin{:}); %Just plot
        end
    else %Nothing to plot so return empty handle
        if nargout==1
            varargout{1}=plot3([],[],[]); %Empty "line array" graphics handle
        end
    end
end

%% Collect output

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
