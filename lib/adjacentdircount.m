function [Lc]=adjacentdircount(L,n)

% function [Lc]=adjacentdircount(L,n)
% ------------------------------------------------------------------------
% This function calculates connectivity (counts the number of neighbouring
% elements) within the input logic L across dimension n. i.e. it that have
% the value 1 in a certain direction. The function can operate on 1D, 2D
% and 3D arrays.
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2008/08/15
% ------------------------------------------------------------------------

%%

L_size=size(L);
if n==3
    n=1;
    L=permute(L,[3,1,2]);
    s=1;
else
    s=0;
end

%%

%Padding array with zeros so that "segmentation" with diff functions
%normally

Lp=zeros(size(L)+2);
if ndims(L)==3
    Lp(2:end-1,2:end-1,2:end-1)=L;
else
    Lp(2:end-1,2:end-1)=L;
end

%Segmentation of groups using diff
dL=diff(Lp,1,n);
IND_starts=find(dL==1); %Starting points where diff is equal to 1
[i_starts,j_starts,k_starts] = ind2sub(size(dL),IND_starts);
IND_ends=find(dL==-1); %Starting points where diff is equal to -1
[i_ends,j_ends,k_ends] = ind2sub(size(dL),IND_ends);

%Getting group lengths (+1 due to diff)
switch n
    case 1
        [j_starts,s_order] = sort(j_starts);
        i_starts=i_starts+1;
        i_starts=i_starts(s_order);
        k_starts=k_starts(s_order);

        [j_ends,s_order] = sort(j_ends);
        i_ends=i_ends(s_order);
        k_ends=k_ends(s_order);

        L_lengths=(i_ends-i_starts)+1;
    case 2
        [i_starts,s_order] = sort(i_starts);
        j_starts=j_starts+1;
        j_starts=j_starts(s_order);
        k_starts=k_starts(s_order);

        [i_ends,s_order] = sort(i_ends);
        j_ends=j_ends(s_order);
        k_ends=k_ends(s_order);

        L_lengths=(j_ends-j_starts)+1;
end

%Creating a matrix of indecises in Lp. Each row corresponds to a group.
%Groups have different lengths but matrix has as many column as the largest
%group. This means that some points appear double.
no_steps=max(L_lengths(:))+(max(L_lengths(:))==1); %Add one if length is 1

[I]=linspacen(i_starts,i_ends,no_steps);
I=floor(I); I=I(~isnan(I));
[J]=linspacen(j_starts,j_ends,no_steps);
J=floor(J); J=J(~isnan(J));
K=linspacen(k_starts,k_ends,no_steps);
K=floor(K); K=K(~isnan(K));

IND = sub2ind(size(Lp),I,J,K);

%Setting up Lc
Lc=zeros(size(Lp));
Lc(IND)=L_lengths*ones(1,no_steps);

%Cropping back to the size of L
if ndims(L)==3
    Lc=Lc(2:end-1,2:end-1,2:end-1); %Cropping back to the size of L
else
    Lc=Lc(2:end-1,2:end-1);
end

if s==1
    Lc=ipermute(Lc,[3,1,2]);
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
