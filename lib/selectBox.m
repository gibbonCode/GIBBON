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