function [Nr,Nr_units]=MRimageUnitScale(N,N_info)

%Scaling data to appropriate units
rescaleSlope=double(N_info.RescaleSlope);
rescaleIntercept=double(N_info.RescaleIntercept);
rescaleType=N_info.RescaleType;

Nr_units=rescaleType;

if isfield(N_info,'ScaleSlope')
    scaleSlope=double(N_info.ScaleSlope);
elseif isfield(N_info,'MRScaleSlope')
    scaleSlope=double(N_info.MRScaleSlope);
else
    scaleSlope=[];
end

if isempty(scaleSlope)
    scaleSlope=1;
    warning('warning, scaleSlope assumed 1!');
end
Nr=((double(N).*rescaleSlope)+rescaleIntercept)./(scaleSlope.*rescaleSlope);



