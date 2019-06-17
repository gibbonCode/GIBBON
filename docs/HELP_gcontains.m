%% gcontains
% Below is a demonstration of the features of the |gcontains| function

%%
clear; close all; clc;

%% Syntax
% |TF=contains(str,strPattern);|

%% Description
% This function is an alternative to the MATLAB function |contains|,
% introduced in MATLAB R2016b. The function attempts to use MATLAB
% |contains| but uses custom implementation if |contains| is not found. See
% also: |contains|.

%% Examples

%% String array

str = ["Mary Ann Jones","Christopher Matthew Burns","John Paul Smith"];

pattern = ["Ann","Paul"];
% TF = contains(str,pattern)
TF = gcontains(str,pattern)

%% Using |IgnoreCase|

str = ["Anne","Elizabeth","Marianne","Tracy"];

pattern = "anne";
% TF = contains(str,pattern,'IgnoreCase',true)

TF = gcontains(str,pattern,'IgnoreCase',true)

%% Character arrays

chr = 'peppers, onions, and mushrooms';

% TF = contains(chr,'onion')
TF = gcontains(chr,'onion')

% TF = contains(chr,'pineapples')
TF = gcontains(chr,'pineapples')