function [cmap_i]=resampleColormap(cmap,n)

% function [cmap_i]=resampleColormap(cmap,n)
% ------------------------------------------------------------------------
%
% This function resamples the input colormap cmap using n steps. Resampling
% is based on linear interpolation. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

ind=(1:1:size(cmap,1))';
ind_i=linspace(1,size(cmap,1),n)';
cmap_i=zeros(n,size(cmap,2));

%Interpolate color data
for q=1:1:size(cmap,2)
    cmap_i(:,q)=interp1(ind,cmap(:,q),ind_i,'linear');
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
