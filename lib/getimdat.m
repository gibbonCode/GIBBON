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
