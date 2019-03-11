function C=makeSymbolicStiffnessTensor(opt)

C=sym(zeros([3 3 3 3]));
I=eye(3,3); 
II1=dyadicProduct(I,I,1);
II3=dyadicProduct(I,I,3);

for i=1:3; 
    for j=1:3;
        for k=1:3; 
            for l=1:3;
                cvar=['c',num2str(i),num2str(j),num2str(k),num2str(l)]; 
         
                
                C(i,j,k,l)=sym(cvar); 
            
            end;
        end; 
    end; 
end;

switch opt
    case 'iso'
        L=(II1+II3)==0;
    case 'transiso'
        
    case'ortho'
        
    case'full'
        L=false(size(C));
    case 'empty'
        L=true(size(C));
end
C(L)=0;

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
