%% imx
% Below is a demonstration of the features of the |imx| function

%%
clear; close all; clc;

%% Syntax
% |hf=imx(varargin);|

%% Description 
% The |imx| function provides a figure window based GUI for 3D image
% segmentation

%% Examples 
%
%% Example: Segmenting MRI data

%% 
% Get a 3D image
load mri;
M=squeeze(D); %example image data set
v=2./[1,1,.4]; %example voxel size

%%
% Start segmentation using |imx|
hf=imx(M,v);

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
