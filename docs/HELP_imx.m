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
% Get a 3D image (see als0 |dcmFolder2MATobject| to import DICOM data)
testCase=5; 
switch testCase
    case 1 %MATLAB brain data
        load mri;
        M=squeeze(D); %example image data set
        v=2./[1,1,.4]; %example voxel size
    case 2 %MRI imported from DICOM files, see also HELP_dcmFolder2MATobject
        defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
        pathName=fullfile(defaultFolder,'data','DICOM','KNEE_UTE');
        loadName=fullfile(pathName,'IMDAT','IMDAT.mat');
        
        IMDAT_struct=load(loadName); %The image data structure 
        G = IMDAT_struct.G; %Geometric/spatial information
        v=G.v; %The voxel size
        M= IMDAT_struct.type_1; %The image data
    case 3
        pathName = '/mnt/data/Dropbox (MIT)/Bryan_Kevin/data/AM/MRI/dcm/Herr_AM_20160321_Herr_AM_20160321/20160321_1.2.276.0.7230010.3.1.2.1635514508.9940.1468594105.272/MR/Left_T1_2';
        loadName=fullfile(pathName,'IMDAT','IMDAT.mat');
        
        IMDAT_struct=load(loadName); %The image data structure 
        G = IMDAT_struct.G; %Geometric/spatial information
        v=G.v; %The voxel size
        M= IMDAT_struct.type_1; %The image data
    case 4
        pathName = '/mnt/data/Dropbox (MIT)/Bryan_Kevin/data/AM/MRI/dcm/Herr_AM_20160321_Herr_AM_20160321/20160321_1.2.276.0.7230010.3.1.2.1635514508.9940.1468594105.272/MR/Right_T1_2';
        loadName=fullfile(pathName,'IMDAT','IMDAT.mat');
        
        IMDAT_struct=load(loadName); %The image data structure
        G = IMDAT_struct.G; %Geometric/spatial information
        v=G.v; %The voxel size
        M= IMDAT_struct.type_1; %The image data
    case 5
        pathName='/mnt/data/Experimental_Data/MRI/AMC/2013/2013_06_14/T1/00801_T1_COR/';
        loadName=fullfile(pathName,'IMDAT','IMDAT.mat');
        
        IMDAT_struct=load(loadName); %The image data structure
        G = IMDAT_struct.G; %Geometric/spatial information
        v=G.v; %The voxel size
        M= IMDAT_struct.type_1; %The image data
end

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