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



 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
