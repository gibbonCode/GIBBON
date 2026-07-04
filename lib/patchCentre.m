function [Vm]=patchCentre(F,V)

% function [Vm]=patchCentre(F,V)
% ------------------------------------------------------------------------
%
%
%
% Change log: 
% 2014
% 2019/04/22 Fixed such that output matches dimentionality of input (output
% used to always be 3D even if input is 2D)
% 2019/04/22 Updated to handle cell input
% ------------------------------------------------------------------------

%%

if isa(F,'cell') %Cell of faces (e.g. potentially of different types)
    Vm=repmat({zeros(size(F,1),3)},size(F)); %Allocate Vm to be a cell like F
    for q=1:1:numel(F) %Loop over face sets       
        Vm{q}=patchCentre(F{q},V); %Override cell entries by face centres
    end    
else
    Vm=zeros(size(F,1),size(V,2)); %Allocate memory for
    for q=1:1:size(V,2)
        X=V(:,q);
        if size(F,1)==1 %Treat single element case since MATLAB is not consistent
            Vm(:,q)=mean(X(F));
        else
            Vm(:,q)=mean(X(F),2);
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
