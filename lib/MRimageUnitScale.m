function [Nr]=MRimageUnitScale(N,N_info)

%Scaling data to appropriate units
rescaleSlope=double(N_info.RescaleSlope);
rescaleIntercept=double(N_info.RescaleIntercept);
% rescaleType=double(N_info.RescaleType);
if isfield(N_info,'ScaleSlope')
    scaleSlope=double(N_info.ScaleSlope);
else
    scaleSlope=[];
end

if isempty(scaleSlope)
    scaleSlope=1;
    warning('warning, scaleSlope assumed 1!');
end
Nr=((double(N).*rescaleSlope)+rescaleIntercept)./(scaleSlope.*rescaleSlope);


