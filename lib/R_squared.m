function [R_sq]=R_squared(yi,fi)  

% function [R_squared]=R_squared()  
% ------------------------------------------------------------------------
% This function calculates R^2 using the data in vector yi and the modelled
% data in the vector fi. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 05/10/2008
% ------------------------------------------------------------------------

y_bar=nanmean(yi);
SS_err=nansum((yi-fi).^2);
SS_tot=nansum((yi-y_bar).^2);
R_sq=1-(SS_err./SS_tot);
 
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
