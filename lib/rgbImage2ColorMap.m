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
