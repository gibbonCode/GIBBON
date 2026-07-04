function [V]=squircle(r,t,n)
% function [V]=squircle(r,t,n)
%-------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------
%%



%Compute x and y components
p=linspace(0,2*pi,n+1)'; p=p(1:end-1); %"angle" parameter

tol=max(eps(p))*2;

if isclose(t,0,tol)
    V=r.*[cos(p) sin(p)];
else    
    %Compute s from tau
    st=r^2*((sqrt(2)/2*(1-t))+t)^2;
    s=(r/st)*sqrt(2*st-r^2);
    
    cos_p=cos(p);
    sin_p=sin(p);
    q=(r/(s*sqrt(2)))*sqrt(1-sqrt(1-s^2*(2*sin_p.*cos_p).^2));
    x=(sign(cos_p))./(abs(sin_p)).*q;
    y=(sign(sin_p))./(abs(cos_p)).*q;
    
    V=real([x y]);
    
    %Treat artifacts
    L=isclose(p,0,tol);
    if any(L)
        V(L,:)=[r 0];
    end
    L=isclose(p,pi/2,tol);
    if any(L)
        V(L,:)=[0 r];
    end
    
    L=isclose(p,pi,tol);
    if any(L)
        V(L,:)=[-r 0];
    end
    
    L=isclose(p,1.5*pi,tol);
    if any(L)
        V(L,:)=[0 -r];
    end
end

%% Flipped quadrants approach

% %Compute x and y components
% if isclose(t,0,tol)
%     p=linspace(0,2*pi,(4*n-4)+1)'; p=p(1:end-1); %angle parameter
%     V=r.*[cos(p) sin(p)];
% else
%     p=linspace(0,pi/2,n)'; p=p(2:end-1); %"angle" parameter
% 
%     %Compute s from tau
%     st=r^2*((sqrt(2)/2*(1-t))+t)^2;
%     s=(r/st)*sqrt(2*st-r^2);
%     
%     cos_p=cos(p);
%     sin_p=sin(p);
%     q=(r/(s*sqrt(2)))*sqrt(1-sqrt(1-s^2*(2*sin_p.*cos_p).^2));
%     x=(sign(cos_p))./(abs(sin_p)).*q;
%     y=(sign(sin_p))./(abs(cos_p)).*q;
%     
%     xx=[r;x;0];
%     yy=[0;y;r];
%    
%     V=real([xx           yy; ...
%            -yy(2:end)    xx(2:end);...
%            -xx(2:end)   -yy(2:end);...
%             yy(2:end-1) -xx(2:end-1)]);        
% end

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
