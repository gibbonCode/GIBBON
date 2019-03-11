function [MU,k]=vonMisesStat(T)

% function [MU,k]=vonMisesStat(T)
% ------------------------------------------------------------------------
% 
% 
% ------------------------------------------------------------------------

%%
MU=angle(mean(exp(1i.*T)));
Rsq=mean(cos(T)).^2+mean(sin(T)).^2;
R=sqrt(Rsq);

p=2; 
ki=R.*(p-Rsq)./(1-Rsq);
diffTol=ki/1000;
qIter=1;
while 1
    kn=ki;
    A=besseli(p/2,ki)./ besseli((p/2)-1,ki);
    ki=ki-((A-R)./(1-A^2-((p-1)/ki)*A));                
    if abs(kn-ki)<=diffTol
        break
    end
    qIter=qIter+1;
end
k=ki;
 
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
