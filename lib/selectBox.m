function [V]=selectBox(varargin)

switch nargin
    case 0
        hf=gcf;    
        zs=0;        
    case 1
        hf=varargin{1};
        zs=0;        
    case 2
        hf=varargin{1};
        zs=varargin{2};        
end

hFunc=get(hf,'WindowButtonMotionFcn');

set(hf,'WindowButtonMotionFcn',{@rbboxFunc,{hf,zs}});

waitforbuttonpress;
    
p1 = get(gca,'CurrentPoint');
p1 = p1(1,1:2);
p1(3)=zs;
p3 = p1;
p2 = p1; p2(1)=p3(1);
p4 = p1; p4(2)=p3(2);
Vr=[p1; p2; p3; p4];

Fr=1:4;

hf.UserData.selectBoxPlot(1)=plotV(Vr([1:end 1],:),'w.-');%patch('Faces',Fr,'Vertices',Vr,'FaceColor','none','EdgeColor',edgeColor,'FaceAlpha',0.1);
hf.UserData.selectBoxPlot(2)=plotV(Vr([1:end 1],:),'k--');

waitforbuttonpress;

p3 = get(gca,'CurrentPoint');
p3 = p3(1,1:2);
p3(3)=zs;

delete(hf.UserData.selectBoxPlot);
hf.UserData=rmfield(hf.UserData,'selectBoxPlot');

V=[p1; p3];

set(hf,'WindowButtonMotionFcn',hFunc);

end

function rbboxFunc(~,~,inputCell)
hf=inputCell{1};
zs=inputCell{2};
hf.UserData.CurrentPoint=get(gca,'CurrentPoint');
p3 = hf.UserData.CurrentPoint;
p3 = p3(1,1:2);
p3(3)=zs;

if isfield(hf.UserData,'selectBoxPlot')
    Vr=[get(hf.UserData.selectBoxPlot(1),'XData')' get(hf.UserData.selectBoxPlot(1),'YData')' get(hf.UserData.selectBoxPlot(1),'ZData')'];
    p1 = Vr(1,:);      
    p2 = p1; p2(1)=p3(1);
    p4 = p1; p4(2)=p3(2);
    Vr=[p1; p2; p3; p4];    
    set(hf.UserData.selectBoxPlot(1),'XData',Vr([1:end 1],1));
    set(hf.UserData.selectBoxPlot(1),'YData',Vr([1:end 1],2));
    set(hf.UserData.selectBoxPlot(1),'ZData',Vr([1:end 1],3));
    set(hf.UserData.selectBoxPlot(2),'XData',Vr([1:end 1],1));
    set(hf.UserData.selectBoxPlot(2),'YData',Vr([1:end 1],2));
    set(hf.UserData.selectBoxPlot(2),'ZData',Vr([1:end 1],3));
%     set(hf.UserData.selectBoxPlot,'Vertices',Vr);
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
