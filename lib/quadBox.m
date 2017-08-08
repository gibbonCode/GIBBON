function [F,V,faceBoundaryMarker]=quadBox(boxDim,boxEl)


%%
F=[]; V=[]; faceBoundaryMarker=[];

[Xtb,Ytb]=meshgrid(0:boxEl(1),0:boxEl(2));

[f,v]=surf2patch(Xtb,Ytb,0*ones(size(Xtb))); %Bottom
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 1*ones(size(f,1),1)];

[f,v]=surf2patch(Xtb,Ytb,boxEl(3)*ones(size(Xtb))); %Top
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 2*ones(size(f,1),1)];

[Yfb,Zfb]=meshgrid(0:boxEl(1),0:boxEl(3));
[f,v]=surf2patch(Yfb,0*ones(size(Yfb)),Zfb); %Front
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 3*ones(size(f,1),1)];

[f,v]=surf2patch(Yfb,boxEl(2)*ones(size(Yfb)),Zfb); %Back
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 4*ones(size(f,1),1)];

[Ylr,Zlr]=meshgrid(0:boxEl(2),0:boxEl(3));
[f,v]=surf2patch(0*ones(size(Ylr)),Ylr,Zlr); %Left
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 5*ones(size(f,1),1)];

[f,v]=surf2patch(boxEl(1)*ones(size(Ylr)),Ylr,Zlr); %Left
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 6*ones(size(f,1),1)];

%%
% Remove double points
[~,ind1,ind2]=unique(V,'rows');
V=V(ind1,:);
F=ind2(F);

%%
% Scale coordinates
maxV=max(V,[],1);
V=V./maxV(ones(size(V,1),1),:);
V=V.*boxDim(ones(size(V,1),1),:);
V=V-boxDim(ones(size(V,1),1),:)/2;
 
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
