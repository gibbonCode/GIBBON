%% spmin
% Below is a demonstration of the features of the |spmin| function

%%
clear; close all; clc;

%% Syntax
% |[minVal,minInd]=spmin(A,B,dim,nanflag,logicRelevant,nanOut);|

%% Description 
% This function is like the min function but is designed for sparse arrays.
% In particular it allows one to "ignore zeros" in the determination of the
% minima. 

%% Examples 
%

%%
% Create example matrix
i=[2 1 1 2  2 3  3 4 4 5  5 5 6 6 7 8];
j=[1 1 2 3  4 5  6 7 8 9 10 11 12 13 13 13];
s=[-1 3 1 2 -1 1 -2 5 5 -1 0 2 3 10 11 NaN];
siz=max([i(:);j(:)]+1)*ones(1,2);
A=sparse(i,j,s,siz(1),siz(2),numel(s));
A=A+A';

full(A) % View matrix

L=sparse(i,j,1,siz(1),siz(2),numel(s));
logicRelevant=(L+L')>0;

%%
% Compute minima allong a certain direction (while omit nan is default)

aMinRows=spmin(A,[],1);
full(aMinRows)

aMinColumns=spmin(A,[],2);
full(aMinColumns)

%%
% Including nans

aMinRows=spmin(A,[],1,'includenan');
full(aMinRows)

aMinColumns=spmin(A,[],2,'includenan');
full(aMinColumns)

%%
% Computing minima across all desired relevant entries (including
% "relevant/real zeros") 

aMinRows=spmin(A,[],1,'omitnan',logicRelevant);
full(aMinRows)

aMinColumns=spmin(A,[],2,'omitnan',logicRelevant);
full(aMinColumns)

%%
% Computin minima across all desired relevant entries and output NaN where
% the sparse array only contains "non-relevant or non-real" zeros. 

nanOut=1;

aMinRows=spmin(A,[],1,'omitnan',logicRelevant,nanOut);
full(aMinRows)

aMinColumns=spmin(A,[],2,'omitnan',logicRelevant,nanOut);
full(aMinColumns)

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
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
