function [D]=getimdat(file_name_IMDAT,opt,d)

D=load(file_name_IMDAT);
IMDAT=D.IMDAT;
siz=IMDAT.dat;
siz_par=IMDAT.par;

switch opt
    case 'dat'
        D=zeros([siz(1) siz(2) siz(3) numel(d)],'uint16');
    case 'par'
        D=repmat(IMDAT.M_info,[siz_par(1) numel(d)]);
end

switch IMDAT.type
    case 'full'
        disp('WARNING! Upload not done. See IMDAT.dat for data');
    case 'split'
        switch opt
            case 'dat'
                for i=1:1:numel(d)
                    IMDAT.load_names{d(i)}
                    load(IMDAT.load_names{d(i)});
                    D(:,:,:,i)=m;
                end
            case 'par'
                for i=1:1:numel(d)
                    load(IMDAT.load_names_par{d(i)});
                    D(:,i)=p;
                end
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
