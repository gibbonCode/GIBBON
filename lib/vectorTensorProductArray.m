function [AB]=vectorTensorProductArray(a,b)

%% 

%Determine multiplication order
if size(a,2)==3 %type B*A
   B=a;
   A=b;
   multType='post';
else %type A*B
    B=b;
    A=a; 
    multType='pre';
end

%Perform product
switch multType
    case 'pre' %type A*B
        %  A1_1*B1 + A1_2*B2 + A1_3*B3
        %  A2_1*B1 + A2_2*B2 + A2_3*B3
        %  A3_1*B1 + A3_2*B2 + A3_3*B3
        AB=[A(:,1).*B(:,1)+A(:,4).*B(:,2)+A(:,7).*B(:,3) ...
            A(:,2).*B(:,1)+A(:,5).*B(:,2)+A(:,8).*B(:,3) ...
            A(:,3).*B(:,1)+A(:,6).*B(:,2)+A(:,9).*B(:,3)];
    case 'post' %type B*A
        %  A1_1*B1 + A1_2*B2 + A1_3*B3
        %  A2_1*B1 + A2_2*B2 + A2_3*B3
        %  A3_1*B1 + A3_2*B2 + A3_3*B3
        AB=[A(:,1).*B(:,1)+A(:,2).*B(:,2)+A(:,3).*B(:,3) ...
            A(:,4).*B(:,1)+A(:,5).*B(:,2)+A(:,6).*B(:,3) ...
            A(:,7).*B(:,1)+A(:,8).*B(:,2)+A(:,9).*B(:,3)];
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
