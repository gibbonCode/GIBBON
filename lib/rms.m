function rms_x=rms(x)

%Calculates RMS for x, NaN's are ignored
x=x(~isnan(x)); %Removing nans
rms_x=sqrt( sum(x.^2)./numel(x) );


end