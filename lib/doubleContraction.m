function C=doubleContraction(varargin)

% function C=doubleContraction(A,B,contractInd)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        A=varargin{1};
        B=varargin{2};
        subContract=[];
    case 3
        A=varargin{1};
        B=varargin{2};
        subContract=varargin{3};
end

if isempty(subContract)
    subContract=[1 2];
end

if numel(subContract)~=2
    error('Two contract indices are needed');
end

subAll=[1 2 3 4];
subNotContract=subAll(~ismember(subAll,subContract));

%% Check input

sizA=size(A);
sizB=size(B);

%Check if A is a fourth-order tensor
siz4=3.*ones(1,4);
siz2=3.*ones(1,2);

%%

if ismatrix(A) && ismatrix(B)    
    C=trace(A*(B.')); %C=Aij.*Bij
elseif ndims(A)==4
    if all(sizA==siz4) && all(sizB==siz2) %4th order and 2nd order tensor
        
        %Creating indices for tensors
        IJKL=gcombvec(1:3,1:3,1:3,1:3)';
        IND_IJKL=sub2indn(siz4,IJKL); %Linear indices for 4th order tensor A
        IND_contract=sub2indn(siz2,IJKL(:,subContract)); %Linear indices for 2nd order tensor B
        
        %Create basis vectors and dyadic product N.B. these are chosen here!
        E=eye(3,3);
        
        %Initialize C
        if isa(A,'double') && isa(B,'double') %if both are doubles
            C=zeros(siz2);
        elseif isa(A,'sym') || isa(B,'sym') %if one of the is a symbolic
            C=sym(zeros(siz2));
        else
            %TO DO! Add warning here
        end
        
        %Sum across indices 
        for q=1:1:numel(IND_IJKL)
            sub1=IJKL(q,subNotContract(1));
            sub2=IJKL(q,subNotContract(2));
            E12=(E(sub1,:).')*E(sub2,:); %dyadic product of basis vectors e1 and e2               
            c=A(IND_IJKL(q)).*B(IND_contract(q)).*E12; %sub-c
            C=C+c; %summed C
        end
    else
        error('A, B have to be two 3x3 arrays or a 3x3x3x3 and 3x3 array')
    end
else
    error('A, B have to be two 3x3 arrays or a 3x3x3x3 and 3x3 array')
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
