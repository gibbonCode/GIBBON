function [E1,E2,E3,G12,G23,G31,v12,v23,v31]=lameInvertHookeOrthotropic(varargin)

%% Parse input

switch nargin
    case 0
        %Initialize all symbolic parameters
        syms mu1 mu2 mu3 lambda11 lambda22 lambda33 lambda12 lambda23 lambda13
    case 9
        %Get parameters symbolic or scalar
        mu1=varargin{1};
        mu2=varargin{2};
        mu3=varargin{3};
        lambda11=varargin{4};
        lambda22=varargin{5};
        lambda33=varargin{6};
        lambda12=varargin{7};
        lambda23=varargin{8};
        lambda13=varargin{9};
end

%%

if any(cellfun(@(x) isa(x,'sym'),varargin))
    symType=1; 
else
    symType=0; 
end

if symType==1
    C_voigt_lame=sym(zeros(6,6));
else
    C_voigt_lame=zeros(6,6); 
end

%Set Lame parameters
C_voigt_lame(1,1)=lambda11+2*mu1;
C_voigt_lame(2,2)=lambda22+2*mu2;
C_voigt_lame(3,3)=lambda33+2*mu3;

C_voigt_lame(1,2)=lambda12;
C_voigt_lame(1,3)=lambda13;
C_voigt_lame(2,1)=lambda12;
C_voigt_lame(3,1)=lambda13;
C_voigt_lame(2,3)=lambda23;
C_voigt_lame(3,2)=lambda23;

C_voigt_lame(4,4)=1/2*(mu2+mu3);
C_voigt_lame(5,5)=1/2*(mu3+mu1);
C_voigt_lame(6,6)=1/2*(mu1+mu2);

%Invert matrix
C_voigt_lame_inv=inv(C_voigt_lame);

%Attempt to simplify if symbolic
if symType==1    
    C_voigt_lame_inv=simplify(C_voigt_lame_inv); %Simplify array
end

%Get alternative parameter set
E1=1/C_voigt_lame_inv(1,1);
E2=1/C_voigt_lame_inv(2,2);
E3=1/C_voigt_lame_inv(3,3);
G12=1/C_voigt_lame_inv(6,6);
G23=1/C_voigt_lame_inv(4,4);
G31=1/C_voigt_lame_inv(5,5);
v12=-(C_voigt_lame_inv(1,2)*E1);
v23=-(C_voigt_lame_inv(2,3)*E2);
v31=-(C_voigt_lame_inv(3,1)*E3);

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
