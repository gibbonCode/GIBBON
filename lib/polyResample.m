function [Vi]=polyResample(V,d_n,interpOpt,interpMethod)

%Parameterize curve using curve length
D = polyCurveLength(V); %Create distance space

if nnz(D>max(eps(V),[],2))<(size(D,1)-1)
    warning('Resampling not possible due to non-unique points (zero distances detected)');
    Vi=[];
else
    
    %Creating even space steps allong curve
    switch interpOpt
        case 1 %d_n is taken to be number of desired points
            D_geo=linspace(0,D(end),d_n);
        case 2 %d_n is taken to be desired point spacing
            d_n_fix=D(end)/round(D(end)/d_n);
            if d_n_fix~=d_n
                disp(['Point spacing was set to: ', num2str(d_n_fix)])
            end
            D_geo=0:d_n_fix:D(end);
    end
    
    %Interpolating curve points
    Vi=zeros(numel(D_geo),size(V,2)); %Initialising Vi
    for qd=1:1:size(V,2)
        Vi(:,qd) = interp1(D,V(:,qd),D_geo,interpMethod);
    end
end




 
%% <-- GIBBON footer text --> 
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
