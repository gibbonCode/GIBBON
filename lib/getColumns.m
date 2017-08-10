function [varargout]=getColumns(V)

% function [varargout]=getColumns(V)
% ------------------------------------------------------------------------
% This function simply outputs each colum in V as a seperate vector. For
% instance if V represents 3 colums of data e.g. X, Y and Z coordinates
% then one could use [X,Y,Z]=getColumns(V);
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/13/08
%------------------------------------------------------------------------
varargout = mat2cell(V,size(V,1),ones(size(V,2),1));


 
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
