function H=errorbarC_XY(x,y,em,ep,w,C,d)

% function H=errorbarC_XY(x,y,em,ep,w,C,d)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%%

n=numel(x);

hold on;
switch d
    case 1
        xu=x;
        xl=x;
        yu=y+ep;
        yl=y-em;

        V=[xu       yu;...
           xl       yl;...
           x-0.5*w  yu;...
           x+0.5*w  yu;...
           x-0.5*w  yl;...
           x+0.5*w  yl];
    case 2
        xu=x+ep;
        xl=x-em;
        yu=y;
        yl=y;

        V=[xu  yu;...
           xl  yl;...
           xu  y-0.5*w;...
           xu  y+0.5*w;...
           xl  y-0.5*w;...
           xl  y+0.5*w];

end

ind1=(n*0+1:n*1)';
ind2=(n*1+1:n*2)';
ind3=(n*2+1:n*3)';
ind4=(n*3+1:n*4)';
ind5=(n*4+1:n*5)';
ind6=(n*5+1:n*6)';

E=[ind1 ind2;...
   ind3 ind4;...
   ind5 ind6];

H=gedge(E,V,repmat(C,3,1));


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
