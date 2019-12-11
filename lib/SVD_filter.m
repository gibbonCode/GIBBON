function [Zm]=SVD_filter(Z,P,T)

%Computer SVD
[U,S,V] = svd(Z);

%Normalising singular values
Ss=diag(S)-min(diag(S));
Ss=Ss./max(Ss);

%Creating smoothening parameters
p_max=P(1); %1=No blurring
p_min=P(2); %0=straight line fit
p=(Ss.*(p_max-p_min))+p_min; %Scale towards singular values

Zm=nan(size(Z));
for i=1:1:size(U,2)
    v=V(:,i); u=U(:,i); s=S(i,i); %components
    if Ss(i)<=T %Filter after threshold
        us = csaps(1:numel(u),u,p(i),1:numel(u))'; %Smooth u
        vs = csaps(1:numel(v),v,p(i),1:numel(v))'; %Smooth v
    else
        vs=v; us=u; %Keep unsmoothened
    end
    z=us*s*vs'; %sub-data
    Zm(:,:,i)=z;
end
Zm=sum(Zm,3);

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
