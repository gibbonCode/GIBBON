%% cellEig
% Below is a demonstration of the features of the |cellEig| function

%% Syntax
% |[V,D]=cellEig(C);|

%% Description 
% Computes eigenvalues and eigenvectors for each matrix contained in the
% cell array C, i.e. [v,d]=eig(c) is executed for each cell entry. The
% output is two cell arrays, i.e. the cell V containing the eigenvectors
% and the cell D containing the eigenvalues. 

%% Examples

%%
close all; clc; clear;

%% Example: Calculating eigenvalues for matrices contained in cells
% Creating example cell containing two matrices

M1=rand(3,3);
M1=M1*M1';
M2=rand(5,5);
M2=M2*M2';

C={M1,M2};
[V,D]=cellEig(C);

%%
% Contained in the output cells are the eigenvectors and eigenvalues of
% each of the matrices e.g. for the first
v1=V{1}
d1=D{1}

%%
% and the second entry
v2=V{2}
d2=D{2}

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
