function [subIndOut,L_valid]=snapSubInd(subIndIn,siz)

%Initialize outputs
L_valid=false(size(subIndIn));
subIndOut=subIndIn;

%Fixing out of range indices
for q=1:1:size(subIndIn,2)
    subInd=subIndIn(:,q); %The current subscript index set
    L_valid(:,q)=(subInd>0) & (subInd<=siz(q)); %Logic for valid indices not out of range
    
    subIndToFix=subInd(~L_valid(:,q)); %Indices to fix
    
    %Snapping out of range indices to boundary
    subIndToFix=(subIndToFix.*(subIndToFix>1))+(subIndToFix<1); %Fix smaller than 1
    subIndToFix=(subIndToFix.*(subIndToFix<=siz(q)))+siz(q).*(subIndToFix>siz(q)); %Fix larger than siz(q)

    %Storing fixed indices in output
    subIndOut_current=subInd;
    subIndOut_current(~L_valid(:,q))=subIndToFix;
    subIndOut(:,q)=subIndOut_current;
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
