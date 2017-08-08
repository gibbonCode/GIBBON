function [A,B,C]=plane_fit(x,y,z)

% function [A,B,C]=plane_fit(x,y,z)
% ------------------------------------------------------------------------
%   Fit a plane to x,y,z data.
%   [A,B,C]=plane_fit(x,y,z) calculates the coefficients A,B,C that fit the data
%   defined by the vectors x,y,z. 
%
%   Uses command svd
%
%   %EXAMPLE: 
%
%         [x,y]=meshgrid(linspace(0,10,20),linspace(0,10,20));
%         a=1; b=2; c=-2;
%         z=(a*x)+(b*y)+c;
%         x=x(:); y=y(:); z=z(:);
%         z=z+(randn(length(z),1));
%         [A,B,C]=plane_fit(x,y,z); 
%         [X,Y]=meshgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%         Z=(A*X)+(B*Y)+C;
%         plot3(x,y,z,'r.'); hold on; grid on;
%         surf(X,Y,Z,'FaceColor','g'); alpha(0.5);
%         title(['a=',num2str(a), ', A=',num2str(A),', b=',num2str(b),', B=',num2str(B),', c=',num2str(c),', C=',num2str(C)]);
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 14/08/2008
% ------------------------------------------------------------------------

P=[mean(x),mean(y),mean(z)];
[U,S,V]=svd([x-P(1),y-P(2),z-P(3)],0);
N=-1/V(end,end)*V(:,end);
A=N(1); B=N(2); C=-P*N;
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
