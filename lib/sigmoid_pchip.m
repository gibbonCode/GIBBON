function [s]=sigmoid_pchip(optStruct)

% function [s]=sigmoid_pchip(optStruct)
%-------------------------------------------------------------------------
% This function generates a sigmoid curve with controllable start and end
% slopes as well as controllable start and end flatness. 
%-------------------------------------------------------------------------

%%

%Default structure
defaultOpt.c1=0;
defaultOpt.c2=defaultOpt.c1;
defaultOpt.r1=0.25;
defaultOpt.r2=defaultOpt.r1;
defaultOpt.n=100;
defaultOpt.nlin=2;
defaultOpt.rMode=1;

[optStruct]=structComplete(optStruct,defaultOpt,0);

c1=optStruct.c1;
c2=optStruct.c2;
r1=optStruct.r1;
r2=optStruct.r2;
nLin=optStruct.nLin;
n=optStruct.n;
rMode=optStruct.rMode;

%%

if rMode==1
    x1=r1;
    x2=r2;
else
    x1=sqrt((r1^2)/(c1^2+1));
    x2=sqrt((r2^2)/(c2^2+1));
end
y1=c1*x1;
y2=(c2*x2);

x2=1-x2;
y2=1-y2;

x1=linspace(0,x1,nLin);
y1=linspace(0,y1,nLin);
x2=linspace(x2,1,nLin);
y2=linspace(y2,1,nLin);

if numel(n)==1
    xi=linspace(0,1,n);
elseif numel(n)>1
    xi=n;
else %n must be empty
    xi=[]; %Make xi empty will force pp form output
end

if ~isempty(xi)
    yi = pchip([x1(:);x2(:)],[y1(:);y2(:)],xi);
else
    pp = pchip([x1(:);x2(:)],[y1(:);y2(:)]);    
end

if numel(n)==1
    s=[xi(:) yi(:)];
elseif numel(n)>1
    s=yi;
else
    s=pp;
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
