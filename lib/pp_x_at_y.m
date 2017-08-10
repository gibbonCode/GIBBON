function x=pp_x_at_y(pp,y,x0,max_iter)

max_fun_eval=max_iter*10;
disp_iter='iter';
OPT_options = optimset( 'MaxFunEvals',max_fun_eval,...
    'MaxIter',max_iter,...
    'TolFun',1e-25,...
    'TolX',1e-25,...
    'Display','off');
[x,OPT_out.resnorm,OPT_out.residual]= fminsearch(@(x) fmin_fvalfind_ppform(pp,x,y),x0,OPT_options);

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
