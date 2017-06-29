function [varargout]=icolorbar(varargin)

hFig=gcf;

switch nargin
    case 0 
        cLim=caxis;
    case 1
        cLim=varargin{1};
end

caxis([cLim(1)-0.5 cLim(2)+0.5]);
hc=colorbar; 
hc.Ticks=cLim(1):1:cLim(2);%linspace(cLim(1)-0.5,cLim(2)+0.5,numel(cLim)+3)
hFig.Colormap=resampleColormap(hFig.Colormap,numel(hc.Ticks));

if nargout==1
    varargout{1}=h;
end