function [cMap]=rgbImage2ColorMap(X,optStruct)

% function [cMap]=rgbImage2ColorMap(X,optStruct)
% ------------------------------------------------------------------------
%
% This function creates a colormap based on an input RGB (red green blue)
% color image (X). Colormap colors are based on the intensities occuring in
% the image. Therefore monotonic color images probably yeild the best
% results. 
% The colormap is harvested from the input image using the settings defined
% in optStruct. The latter contains the fields: 
% n: The number of colormap levels (default 250)
% normFactor: The normalisation factor e.g. 255 (default is maximum value
% occuring in image).
% numBins: Number of bins to sample the image with (colors are averaged for
% each bin, default=n)
%
% See also: |colormap|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/04/01 Added to GIBBON
%------------------------------------------------------------------------

%% Parse input

%Get options from structure or use defaults

if isfield(optStruct,'n')
   n=optStruct.n; 
else
    n=250;
end

if isfield(optStruct,'normFactor')
   normFactor=optStruct.normFactor; 
else
    normFactor=max(X(:));
end

if isfield(optStruct,'numBins')
   numBins=optStruct.numBins; 
else
    numBins=n;
end

if isfield(optStruct,'imageIntensityLabel')
   imageIntensityLabel=optStruct.imageIntensityLabel; 
else
    imageIntensityLabel=mean(X,3);
end

%Convert image to a double if it isn't already
if ~isa(X,'double')
    X=double(X);
end

%% Sample colors for desired number of intensity levels

%Compute normalised intensity in numBins+1 integers
imageIntensityLabel=imageIntensityLabel-min(imageIntensityLabel(:)); 
imageIntensityLabel=imageIntensityLabel./max(imageIntensityLabel(:));
imageIntensityLabel=imageIntensityLabel.*numBins; 
imageIntensityLabel=round(imageIntensityLabel); 

intAll=unique(imageIntensityLabel(:))';

%Get color layers
R=X(:,:,1);
G=X(:,:,2);
B=X(:,:,3);

cMap_sub=zeros(numel(intAll),3);
for q=1:1:numel(intAll)
    indNow=find(imageIntensityLabel==intAll(q));
    cMap_sub(q,:)=[mean(R(indNow)) mean(G(indNow)) mean(B(indNow))];
end

%% Interpolate to evenly across intensities with n levels

intRange=linspace(0,numBins,n)';

cMap=zeros(n,3);
for q=1:1:3;
    cMap(:,q)=interp1(intAll,cMap_sub(:,q),intRange,'linear');
end

%% Normalise colors
cMap=cMap./normFactor; 

%%
if any(cMap(:)>1) || any(cMap(:))<0
    error('Colormapped values should be in the range [0 1]. Alter input image and/or normFactor');
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
