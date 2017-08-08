function [varargout]=icolorbar(varargin)

hFig=gcf;

switch nargin
    case 0 
        cLim=caxis;
    case 1
        cLim=varargin{1};
end

caxis([cLim(1)-0.5 cLim(2)+0.5]);
hc=colorbar; 
hc.Ticks=cLim(1):1:cLim(2);%linspace(cLim(1)-0.5,cLim(2)+0.5,numel(cLim)+3)
hFig.Colormap=resampleColormap(hFig.Colormap,numel(hc.Ticks));

if nargout==1
    varargout{1}=h;
end
 
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
