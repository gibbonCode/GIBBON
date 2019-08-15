function hf=sv2(varargin)

%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
        viewOpt=[];
    case 2
        M=varargin{1};
        v=varargin{2}; 
        viewOpt=[];
    case 3
        M=varargin{1};
        v=varargin{2};
        viewOpt=varargin{3};
end

%%

M=double(M);

if ~isempty(viewOpt)
    cLim=viewOpt.cLim;    
    logicMask=viewOpt.logicMask;     
else
    cLim=[min(M(:)) max(M(:))];
    logicMask=true(size(M)); 
end
cLimScaled=round((cLim./max(M(:)))*100);

M(~logicMask)=NaN;

%%
% Plot settings
fontSize=10;
fontColor='w';
cMap=gray(250);

% figStruct.vcw=0; %Currently not compatible with vcw
figStruct.Color='k'; %Figure background color
figStruct.ColorDef='black'; %Setting colordefinitions to black
% figStruct.ScreenOffset=20; %Setting spacing of figure with respect to screen edges

%%

%Defining row, column and slice indicices for slice patching
sliceIndexI=round(size(M,1)/2); %(close to) middle row
sliceIndexJ=round(size(M,2)/2); %(close to) middle column
sliceIndexK=round(size(M,3)/2); %(close to) middle slice

%%

[ax,ay,az]=im2cart([size(M,1)+1 0],[size(M,2)+1 0],[size(M,3)+1 0],v);
axLim=[ax(2) ax(1) ay(2) ay(1) az(2) az(1)];

hf=cFigure(figStruct);

for q=1:1:4
    subplot(2,2,q);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    hold on;
    axis equal; axis tight; axis vis3d; axis(axLim); grid on; box on;
    
    switch q
        case 1
            title('JK-view','color','r');
            view(0,0);
        case 2
            title('IK-view','color','g');
            view(-90,0);
        case 3
            title('IJ-view','color','b');
            view(0,90);
        case 4
            title('IJK-view ','color',fontColor);
            colorbar;
            view(3);
    end
    colormap(cMap);
    caxis(cLim);
    set(gca,'fontSize',fontSize);
    H(q)=gca;
end
drawnow;

%% Set user data

hf.UserData.sv2.M=M;
hf.UserData.sv2.siz=size(M);
hf.UserData.sv2.v=v;
hf.UserData.sv2.patchTypes={'si','sj','sk'};
hf.UserData.sv2.H=H;
hf.UserData.sv2.hp=NaN(1,6);
hf.UserData.sv2.hpp=NaN(1,6);
hf.UserData.sv2.sliceIndices=[sliceIndexI sliceIndexJ sliceIndexK];
hf.UserData.sv2.axLim=axLim;
hf.UserData.sv2.cLim=cLim;
hf.UserData.sv2.fontColor=fontColor;

updateSlices(hf);
 
%%

set(hf,'WindowButtonDownFcn', {@ButtonDownFunction,hf},'BusyAction','cancel');

%%

drawnow;
end

function ButtonDownFunction(~,~,hf)
v=hf.UserData.sv2.v;

cax = overobj2('axes');

if isempty(cax)
    %     cax=gca; %this gets current axis or if none exists creates one
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end

if isempty(cax)
    return
end

%Key actions
pt_ax = get(cax, 'CurrentPoint');

if cax==hf.UserData.sv2.H(1) %JK view -> I slice
    %     hp=plot(pt_ax(1,1),pt_ax(1,3),'r.');
    [~,jj,kk]=cart2im(pt_ax(1,1),0,pt_ax(1,3),v);
    jj=round(jj);
    kk=round(kk);
    hf.UserData.sv2.sliceIndices(2)=jj;
    hf.UserData.sv2.sliceIndices(3)=kk;
    hf=fixSliceIndices(hf);
elseif cax==hf.UserData.sv2.H(2) %IK view -> J slice
    %     hp=plot(pt_ax(1,1),pt_ax(1,3),'r.');
    [ii,~,kk]=cart2im(0,pt_ax(1,2),pt_ax(1,3),v);
    ii=round(ii);
    kk=round(kk);
    hf.UserData.sv2.sliceIndices(1)=ii;
    hf.UserData.sv2.sliceIndices(3)=kk;
    hf=fixSliceIndices(hf);
elseif cax==hf.UserData.sv2.H(3) %IJ view -> K slice
    %     hp=plot(pt_ax(1,1),pt_ax(1,2),'r.');
    [ii,jj,~]=cart2im(pt_ax(1,1),pt_ax(1,2),0,v);
    ii=round(ii);
    jj=round(jj);
    hf.UserData.sv2.sliceIndices(1)=ii;
    hf.UserData.sv2.sliceIndices(2)=jj;
    hf=fixSliceIndices(hf);
else     
    return
end

updateSlices(hf);

end

function hf=fixSliceIndices(hf)

siz=hf.UserData.sv2.siz; %Get image size
sliceIndices= hf.UserData.sv2.sliceIndices; %Get proposed slice indices

%Fix indices to be between 1 and siz
sliceIndices(sliceIndices<1)=1; 
sliceIndices(sliceIndices>siz)=siz(sliceIndices>siz);
hf.UserData.sv2.sliceIndices=sliceIndices; 

end

function updateSlices(hf)

sliceIndices = hf.UserData.sv2.sliceIndices;
M=hf.UserData.sv2.M;
v=hf.UserData.sv2.v;
axLim=hf.UserData.sv2.axLim;
[xx,yy,zz]=im2cart(sliceIndices(1),sliceIndices(2),sliceIndices(3),v);

for dirOpt=1:1:3
    patchType=hf.UserData.sv2.patchTypes{dirOpt};
    logicPatch=false(size(M));
    sliceIndex=sliceIndices(dirOpt);

    switch dirOpt
        case 1
            logicPatch(sliceIndex,:,:)=1;
        case 2
            logicPatch(:,sliceIndex,:)=1;
        case 3
            logicPatch(:,:,sliceIndex)=1;
    end
        
    if isnan(hf.UserData.sv2.hp(dirOpt))        
        [F,V,C]=im2patch(M,logicPatch,patchType,v);    
        subplot(2,2,dirOpt);
        hf.UserData.sv2.hp(dirOpt)= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','none');
        
        subplot(2,2,4);
        hf.UserData.sv2.hp(dirOpt+3)= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','none');
    else
        set(hf.UserData.sv2.hp(dirOpt),'CData',M(logicPatch));
        set(hf.UserData.sv2.hp(dirOpt+3),'CData',M(logicPatch));
        V=get(hf.UserData.sv2.hp(dirOpt),'Vertices');
        switch dirOpt
            case 1
                V(:,2)=(sliceIndex-0.5).*v(1);
            case 2
                V(:,1)=(sliceIndex-0.5).*v(2);
            case 3
                V(:,3)=(sliceIndex-0.5).*v(3);
        end
        set(hf.UserData.sv2.hp(dirOpt),'Vertices',V);        
        set(hf.UserData.sv2.hp(dirOpt+3),'Vertices',V);        
    end
    
    subplot(2,2,dirOpt);    
    
    if any(isnan(hf.UserData.sv2.hpp))        
        f=1:4;
        switch dirOpt
            case 1
                vv=[xx axLim(3) axLim(6); xx axLim(4) axLim(6); xx axLim(4) axLim(5); xx axLim(3) axLim(5);];
                hf.UserData.sv2.hpp(1)=patch('faces',f,'vertices',vv,'faceColor','g','EdgeColor','g','faceAlpha',0.1);
                vv=[axLim(1) axLim(3) zz; axLim(2) axLim(3) zz; axLim(2) axLim(4) zz; axLim(1) axLim(4) zz;];
                hf.UserData.sv2.hpp(2)=patch('faces',[1 2 3 4],'vertices',vv,'faceColor','b','EdgeColor','b','faceAlpha',0.1);
            case 2
                vv=[axLim(1) yy axLim(5); axLim(2) yy axLim(5); axLim(2) yy axLim(6); axLim(1) yy axLim(6);];
                hf.UserData.sv2.hpp(3)=patch('faces',f,'vertices',vv,'faceColor','r','EdgeColor','r','faceAlpha',0.1);
                vv=[axLim(1) axLim(3) zz; axLim(2) axLim(3) zz; axLim(2) axLim(4) zz; axLim(1) axLim(4) zz;];
                hf.UserData.sv2.hpp(4)=patch('faces',f,'vertices',vv,'faceColor','b','EdgeColor','b','faceAlpha',0.1);
            case 3
                vv=[axLim(1) yy axLim(5); axLim(2) yy axLim(5); axLim(2) yy axLim(6); axLim(1) yy axLim(6);];
                hf.UserData.sv2.hpp(5)=patch('faces',f,'vertices',vv,'faceColor','r','EdgeColor','r','faceAlpha',0.1);
                vv=[xx axLim(3) axLim(6); xx axLim(4) axLim(6); xx axLim(4) axLim(5); xx axLim(3) axLim(5);];
                hf.UserData.sv2.hpp(6)=patch('faces',f,'vertices',vv,'faceColor','g','EdgeColor','g','faceAlpha',0.1);
        end
    else        
        switch dirOpt
            case 1                
                set(hf.UserData.sv2.hpp(1),'Vertices',[xx axLim(3) axLim(6); xx axLim(4) axLim(6); xx axLim(4) axLim(5); xx axLim(3) axLim(5);]);                                
                set(hf.UserData.sv2.hpp(2),'Vertices',[axLim(1) axLim(3) zz; axLim(2) axLim(3) zz; axLim(2) axLim(4) zz; axLim(1) axLim(4) zz;]);                                
            case 2
                set(hf.UserData.sv2.hpp(3),'Vertices',[axLim(1) yy axLim(5); axLim(2) yy axLim(5); axLim(2) yy axLim(6); axLim(1) yy axLim(6);]);
                set(hf.UserData.sv2.hpp(4),'Vertices',[axLim(1) axLim(3) zz; axLim(2) axLim(3) zz; axLim(2) axLim(4) zz; axLim(1) axLim(4) zz;]);
            case 3
                set(hf.UserData.sv2.hpp(5),'Vertices',[axLim(1) yy axLim(5); axLim(2) yy axLim(5); axLim(2) yy axLim(6); axLim(1) yy axLim(6);]);
                set(hf.UserData.sv2.hpp(6),'Vertices',[xx axLim(3) axLim(6); xx axLim(4) axLim(6); xx axLim(4) axLim(5); xx axLim(3) axLim(5);]);
        end        
    end
end

titleString=['Image coordinates: ',num2str(sliceIndices(1)),' ',num2str(sliceIndices(2)),' ',num2str(sliceIndices(3)),', Color limits: ',sprintf('%.2f %.2f',hf.UserData.sv2.cLim)];
hf.Name=titleString;
drawnow;

end

function setColorLimits(~,~,inputCell)

hf=inputCell{1};
jSlider_C=inputCell{2};

cLimScaled(1) = get(jSlider_C,'LowValue');
cLimScaled(2) = get(jSlider_C,'HighValue');

hf.UserData.sv2.cLim=(cLimScaled/100).*max(hf.UserData.sv2.M(:));

for q=1:1:numel(hf.UserData.sv2.H)
    hf.UserData.sv2.H(q).CLim=hf.UserData.sv2.cLim;
end

sliceIndices = hf.UserData.sv2.sliceIndices;
titleString=['Image coordinates: ',num2str(sliceIndices(1)),' ',num2str(sliceIndices(2)),' ',num2str(sliceIndices(3)),', Color limits: ',sprintf('%.2f %.2f',hf.UserData.sv2.cLim)];
hf.Name=titleString;

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
