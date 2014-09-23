function [F]=fbg2force(FBG_temp,FBG_temp_zero,FBG_strain,FBG_strain_zero)


%% CALCULATING FBG FIBRE STRAIN

Se=0.7773; %7.773e-7 1/microstrain, or 0.7773/strain
FBG_e=(1/Se).*(log(FBG_strain/FBG_strain_zero)-log(FBG_temp/FBG_temp_zero));

%% GETTING CALIBRATION DATA

load_path='C:\Users\kmmoerman\00_WORK\04_EXPERIMENTAL_DATA\04_FBG\CALIBRATION\2009_11_26_calibration_01\CALIBRATION_RESULTS\';
file_name=['CALI_FIT_ALL.mat'];
load([load_path,file_name]);

%% CALCULATING FORCE

F=ppval(CALI_PP_ALL,FBG_e);

end
