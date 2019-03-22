function hg=gauss_kernel(k,nd,f,m)

% function hg=gauss_kernel(k,nd,f,m)
% ------------------------------------------------------------------------
% 
%
% ------------------------------------------------------------------------

%%
switch m
    case {'sigma',1}
        %NOTE: This method is equivalent to using hg = fspecial('gaussian',[k k], S)
        S=f3;
        switch nd
            case 1
                x = linspace(-((k-1)/2),((k-1)/2),k);
                hg=exp(-(x.^2)./(2*S^2));
            case 2
                [x,y] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
                hg=exp(-(x.^2 + y.^2)./(2*S^2));
            case 3
                [x,y,z] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
                hg=exp(-(x.^2 + y.^2 + z.^2)./(2*S^2));
        end
    case {'width',2}
        % Here "f" defines "where the bell curve is" at the edges of the
        % kernel e.g. f=2 results in twice the standard deviation.
        switch nd
            case 1
                x = linspace(-f,f,k);
                hg=exp(-(x.^2)./2);
            case 2
                [x,y] = meshgrid(linspace(-f,f,k));
                hg=exp(-(x.^2 + y.^2)./2);
            case 3
                [x,y,z] = meshgrid(linspace(-f,f,k));
                hg=exp(-(x.^2 + y.^2 + z.^2)./2);
        end
    otherwise
        warning('False input for method! Valid options are "sigma" and "width"');
end

hg(hg<eps*max(hg(:)))=0;

%Making sure mask sums to 1
hg=hg./sum(hg(:));


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
