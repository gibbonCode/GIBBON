function Fmin=fmin_symbolic_form(Pn,OPT_struct,F)

% %Getting function variables 
% for i=1:1:numel(P)    
%     eval(['P',num2str(i),'=P(i);']);
% end

X=OPT_struct.FX;

%Evaluate function
P=sym('P',size(Pn)); %Create P as symbolic
Ff=subs(subs(F,P,Pn)); %Substitution of P and other

%Derive optimisation function value
switch OPT_struct.method
    case 'fminsearch'
        Fmin=sum((Ff(:)-OPT_struct.FY(:)).^2);
    case 'lsqnonlin' 
        Fmin=(Ff-OPT_struct.FY);
end
Fmin=double(Fmin);

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
