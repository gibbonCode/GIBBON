function maximizeFigureWindow(varargin)

%% Parse input

switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

if isempty(hf)
    hf=gcf; %Take the current figure window handle if the handle is empty
end

%% Turn off warning (JavaFrame will become obsolete, this will lead to error in the future)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

%% Maximize using JaveFrame
try
    pause(0.01);     
    jFrame=get(hf,'JavaFrame');
%     set(jFrame,'Maximized',1);
    jFrame.setMaximized(1);
catch exception    
    rethrow(exception)
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
