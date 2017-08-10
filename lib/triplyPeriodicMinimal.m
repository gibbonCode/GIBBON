function S=triplyPeriodicMinimal(varargin)

%Parse input
switch nargin
    case 2
        P=varargin{1};
        X=P(:,1); Y=P(:,2); Z=P(:,3); %Get coordinates
        typeStr=varargin{2};
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        typeStr=varargin{4};
    otherwise
        error('False number of input arguments.');
end

%Evaluate metric on coordinates
switch typeStr
    case 'p' %Schwarz P
        S=(cos(X)+cos(Y)+cos(Z));
    case 'd' %Schwarz D
        S=(sin(X).*sin(Y).*sin(Z))+(sin(X).*cos(Y).*cos(Z))+(cos(X).*sin(Y).*cos(Z))+(cos(X).*cos(Y).*sin(Z));
    case 'g' %Schoen Gyroid
        S=(sin(X).*cos(Y))+(sin(Y).*cos(Z))+(cos(X).*sin(Z)); %Schoen Gyroid
    case 'n' %Neovius
        S=3*(cos(X)+ cos(Y)+ cos(Z))+ (4*cos(X).*cos(Y).*cos(Z));
    case 'w'
        S=cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X)-cos(X).*cos(Y).*cos(Z);
    case 'pw'
        S=(4.*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-3.*cos(X).*cos(Y).*cos(Z))+2.4;
    otherwise
        error('unknown surface type requested')
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
