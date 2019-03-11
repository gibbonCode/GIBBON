function A=numReplace(A,a,b)

% function A=numReplace(A,a,b)
% ------------------------------------------------------------------------
% Replacing the entries in a occuring in A with the entries in b. 
% 
% Example: 
% %Define input array
% A=[0,-6,3,0,0;...
%     -1,-4,-7,-9,-9;...
%     -4,-7,11,-5,12;...
%     10,-7,5,-7,13;...
%     0,11,-2,-6,2;...
%     12,4,2,-5,NaN];
% a=[2 -5 nan 0]; %Entries to replace
% b=[991 992 993 994]; %Entries to take their place
% B=numReplace(A,a,b); % Replacing the numbers using |numReplace|
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/01
%------------------------------------------------------------------------

if numel(a)~=numel(b)
    error('numel(a)~=numel(b)')
end

if numel(a)~=numel(unique(a))
    error('numel(a)~=numel(unique(a))')
end

%Dealing with NaN entries
L_nan=isnan(A); %Logic for the nan entries
nanRep=max(A(~L_nan))+1; %Value to replace nan's with (arbitrary number not in A)
A(L_nan)=nanRep; 
a(isnan(a))=nanRep;

%Sorting input vectors;
[a,I]=sort(a); %Values to replace
b=b(I); %Numbers to appear in output in place of the entries in a

%Replacing numbers in a with those in b
L_mem=ismember(A,a); %Finding members common to A and a
valRep=A(L_mem); %These are the (non-unique) entries to be replaces
[~,~,indReplace]=unique(valRep); %Use to unique operation to get indices into b
A(L_mem)=b(indReplace); %Now replace the entries
 
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
